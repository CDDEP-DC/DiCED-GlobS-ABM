using LinearAlgebra
include("utils.jl")
include("fileutils.jl")

## remote workers need these definitions
##
## NOTE : remote workers currently don't have utils.jl (to save memory, loading libraries)
##
@everywhere begin

using Random
using SparseArrays
using Distributions

## samples with replacement from collection a
## if a is empty, returns an empty vector instead of an error
function rsamp(a::AbstractArray{T}, n::I)::Vector{T} where {T<:Any,I<:Integer}
    return isempty(a) ? T[] : rand(a,n)
end

##
## can create any type of event, just need to write a handler for it like event handlers below
## can create any agent states, just need to create events+handlers to transition between states
##
## note, we don't actually want these structs to be mutable, but immutable structs with equal
##   fields are "identical", and PriorityQueue doesn't like that
##  (using immutable structs with an autoincremented field to distinguish them is worse)
mutable struct infectionEvent <: simEvent
    const t::UInt32
    const agentid::UInt32
end

mutable struct becomeContagious <: simEvent
    const t::UInt32
    const agentid::UInt32
end

mutable struct becomeRecovered <: simEvent
    const t::UInt32
    const agentid::UInt32
end

## periodStuff handles out-of-network infections (people who work or live outside of the synth area)
mutable struct periodicStuff <: simEvent
    const t::UInt32
end

## event to trigger reporting
mutable struct reportingEvent <: simEvent
    const t::UInt32
end

## separate handler for infecting people while on holiday
mutable struct becomeHolidayContagious <: simEvent
    const t::UInt32
    const agentid::UInt32
end

## struct that stores data within a SimUnit
## (this could have just been a single dict, but this way the types are identified, helps the compiler)
mutable struct SimData <: abstractSimData
    netw::Dict{Symbol,SparseMatrixCSC{Bool,UInt32}} ## together in a dict so they can be treated as a collection
    sets::Dict{Symbol,Set{UInt32}} ## which agents belong to this sim unit, along with any subsets
    inf_probs::Dict{Symbol,Float64} ## store in a dict so a different one can be associated with each network
    distr_fns::Dict{Symbol,Function} ## ditto for fns used to determine contact events
    distr_params::Dict{Symbol,Tuple} ## params for those fns
    loc_lookup::Dict{Symbol,Dict{UInt32,UInt32}} ## for indexing ephemeral contact networks by location

    n_agents::Int
    I_set::Set{UInt32} ## current infected
    R_set::Set{UInt32} ## current recovered
    t_recovery::UnitRange{Int64} ## range of recovery times
    report_freq::UInt32
    periodic_stuff_period::UInt32 ## mainly for updating out-of-network infections

    dummies_assigned::Vector{UInt32} ## workers without households; special logic needed; vector because randsubseq can't handle Set
    outw_assigned::Vector{UInt32} ## workers without workplaces
    mean_hh_connections::Float64 ## mainly for the dummies
    mean_wp_connections::Float64 ## mainly for out-workers

    report_groups::Vector{Symbol} ## sets for which data should be reported
    cumI::Dict{Symbol,Int} ## cumulative infection counts for each set
    report_series::Dict{Symbol,Vector{Int}} ## dict of time series for reporting
    geo_assigned::Dict{UInt32,UInt16} ## for tracking # infections by home location
    cumI_by_geo::Vector{Int} ## each location is a vector index
    report_by_geo::Vector{Vector{Int}} ## time series by geo; each location is a vector index

    intervals::Dict{Symbol, UnitRange{Int64}} ## time intervals for special treatment; holidays etc
    flags::Set{Symbol} ## misc boolean flags, mainly for testing
    ## constructor returns an empty struct; filled using external inputs by init_sim_unit fn below
    SimData() = new()
end

##
## event handlers, dispatch on event type
## worker loop expects these functions to return future events as a Vector{simEvent}
##

## most efficient way to do logging/reporting:
## add code to existing events (e.g., increment case count in E->I event)
##  and queue a "data summary" event similar to sync event

function handle_event!(u::SimUnit, e::infectionEvent)::Vector{simEvent}
    ## this event could have occurred in another simunit, in the past; that's ok because nothing has happened yet
    ## just return an event for the target to become contagious at the correct time
    ## (note if an agent becomes infected twice, the earlier one will take effect & the later will do nothing)
    effect_time::UInt32 = e.t + u.t_inc
    #if in(effect_time, u.d.intervals[:holiday])
    #    return [becomeHolidayContagious(effect_time, e.agentid)]
    #else
    return [becomeContagious(effect_time, e.agentid)]
    #end
end

##
## functions used in becomeContagious handler to specify the distribution of contacts
##

## all neighbors are contacted once or with equal frequency/intensity
function distr_all(neigh::AbstractVector{I}, params::T)::Vector{UInt32} where {I<:Integer, T<:Tuple}
    return neigh
end

## N random contact events; samples with replacement from neigh; N is the first/only value in params
##   note, params must be a single-Int tuple (N,) if using this fn 
function distr_const(neigh::AbstractVector{I}, params::Tuple{Int64})::Vector{UInt32} where I<:Integer
    return rsamp(neigh, params[1])
end

## Poission(L) gives the number of contact events, then sample with replacement; L is the first/only value in params
function distr_pois(neigh::AbstractVector{I}, params::Tuple{Int64})::Vector{UInt32} where I<:Integer
    n = rand(Poisson(params[1]))
    return rsamp(neigh, n)
end

## Geometric (discrete analogue of Exponential) gives the number of contact events; params[1] is the mean
function safe_rgeo(u::R) where R<:Real
    return u > 0 ? rand(Geometric(1/(1+u))) : 0 ## must be > 0
end
function distr_exp(neigh::AbstractVector{I}, params::Tuple{Int64})::Vector{UInt32} where I<:Integer
    n = safe_rgeo(params[1])
    return rsamp(neigh, n)
end

## fn to indicate no contacts of a certain type
function distr_zero(neigh::AbstractVector{I}, params::T)::Vector{UInt32} where {I<:Integer, T<:Tuple}
    return UInt32[]
end

## if keeping overall contact rate constant, don't allow mean loc contacts to be greater than # neighbors
function w_test_params(d::SimData,i::Integer)
    N_work_neigh = nnz(view(d.netw[:wp],:,i))
    params_loc_work = (get(d.distr_params[:loc_work],1,0) > N_work_neigh) ? 
        (N_work_neigh, d.distr_params[:loc_work][2:end]...) : 
        d.distr_params[:loc_work]
    params_nonloc = (get(d.distr_params[:nonloc],1,0) > N_work_neigh) ?
        (N_work_neigh, d.distr_params[:nonloc][2:end]...) :
        d.distr_params[:nonloc]
    return (params_loc_work, params_nonloc)
end
function h_test_params(d::SimData,i::Integer)
    N_res_neigh = nnz(view(d.netw[:hh],:,i)) + nnz(view(d.netw[:gq],:,i))
    params_loc_res = (get(d.distr_params[:loc_res],1,0) > N_res_neigh) ?
        (N_res_neigh, d.distr_params[:loc_res][2:end]...) :
        d.distr_params[:loc_res]
    params_nonloc = (get(d.distr_params[:nonloc],1,0) > N_res_neigh) ?
        (N_res_neigh, d.distr_params[:nonloc][2:end]...) :
        d.distr_params[:nonloc]
    return (params_loc_res,params_nonloc)
end

## if keeping overall contact rate constant, don't allow mean loc contacts to be greater than # neighbors
function get_dparams(d::SimData, i::UInt32, keynames::Vector{Symbol})::Vector{Tuple}
    tmp_copy = Dict{Symbol,Tuple}(k=>d.distr_params[k] for k in keynames)
    if in(:w_test, d.flags)
        tmp_copy[:loc_work], tmp_copy[:nonloc] = w_test_params(d,i)
    end
    if in(:h_test, d.flags)
        tmp_copy[:loc_res], tmp_copy[:nonloc] = h_test_params(d,i)
    end
    return [tmp_copy[k] for k in keynames]
end

function resists_infection(i::UInt32, e::simEvent, d::SimData)::Bool
    return in(i, d.I_set) || in(i, d.R_set) ## testing set membership is O(1)
end

function loc_res(d::SimData, i::UInt32)::UInt32
    return get(d.loc_lookup[:res], i, UInt32(0))
end
function loc_work(d::SimData, i::UInt32)::UInt32
    return get(d.loc_lookup[:work], i, UInt32(0))
end

## TODO
function loc_holiday(d::SimData, i::UInt32)::UInt32
    return UInt32(0)
end

function update_counts!(d::SimData, i::UInt32)
    for k in d.report_groups
        if i in d.sets[k] ## agent id belongs to this group
            d.cumI[k] += 1
        end
    end
    ## update geo vector
    g = d.geo_assigned[i]
    if g > 0 ## geo_assigned is 0 for dummies, who don't have a home location
        d.cumI_by_geo[g] += 1
    end
end

function network_infections(netw::A, col_idx::UInt32, distr_fn::F, distr_params::T, p_inf::Float64)::Vector{UInt32} where {A<:AbstractArray, F<:Function, T<:Tuple}
    if col_idx > 0 ## index exists
        if distr_fn == distr_zero ## avoid doing unnecessary work
            return UInt32[]
        else
            neigh = findall(view(netw,:,col_idx)) ## note, findall() on a sparsearray is not really a "find", it's an O(1) lookup
            contacts = distr_fn(neigh, distr_params)
            return randsubseq(contacts,p_inf) ## randsubseq promises efficient Bernouilli sampling
        end
    else
        return UInt32[]
    end
end

##
## the becomeContagious handler
##   if agent is susceptible, change state, and generate all future infection events that result
##
function handle_event!(u::SimUnit, e::becomeContagious)::Vector{simEvent}
    d = u.d
    i = e.agentid
    if resists_infection(i, e, u.d)
        return simEvent[] ## agent not susceptible
    else
        push!(d.I_set, i) ## append to current infected set
        update_counts!(u.d, i)
        duration = rand(d.t_recovery)
        recov_event = becomeRecovered(e.t + duration, i)
        ## collect references to networks defined in the simunit (we will broadcast the infection fn over these)
        networks = [d.netw[k] for k in [:hh, :wp, :sch, :gq, :loc_matrix_res, :loc_matrix_work]]
        ## column = agent id, except location networks where it's a location index, and nonlocal "network" that's just 1 column
        col_idxs = [i, i, i, i, loc_res(u.d,i), loc_work(u.d,i)]
        ## distribution function for # of contacts in each network; currently using non-hh for gq's
        distr_fns = [d.distr_fns[k] for k in [:hh, :non_hh, :non_hh, :non_hh, :loc_res, :loc_work]]
        distr_params = get_dparams(u.d,i,[:hh, :non_hh, :non_hh, :non_hh, :loc_res, :loc_work])        
        inf_ps::Vector{Float64} = [d.inf_probs[k] for k in [:p_inf_hh, :p_inf, :p_inf, :p_inf, :p_inf_loc, :p_inf_loc]]

        ## make modifications based on whichever special time intervals are currently active
        active_keys::Vector{Symbol} = [k for (k, i_range) in d.intervals if in(e.t, i_range)]
        for i_key in active_keys
            ## switch on i_key
            if (i_key == :nonessential_wp_closed)
                ## workplaces are closed except low income; disable wp and wp loc networks, unless agent is in low inc wp
                in(i, d.sets[:low_inc_workplace]) || (col_idxs[[2,6]] .= 0) ## set col idx to 0 to exclude the network
            
            elseif (i_key == :sch_closed)
                ## schools are closed; disable school network for students and work networks for teachers
                col_idxs[3] = 0 ## school network only exists for k12 students
                in(i, d.sets[:k12_worker]) && (col_idxs[[2,6]] .= 0)
            end
        end

        ## f.() syntax broadcasts a fn over collections; returns a vector of f's return vals; those are vectors, so join them with reduce(vcat)
        ##   unique() because infecting someone twice doesn't have any additional effect
        infected::Vector{UInt32} = unique(reduce(vcat, network_infections.(networks, col_idxs, distr_fns, distr_params, inf_ps)))
        infection_events = [infectionEvent(e.t + rand(0:duration), targ) for targ in infected]
        ##
        ## TODO: log number of secondary infections
        ##

        ## return infection events and recovery event
        return [infection_events; recov_event]
    end
end


## become contagious while on holiday
##  cheat to keep the overall infection rate unchanged
function handle_event!(u::SimUnit, e::becomeHolidayContagious)::Vector{simEvent}
    d = u.d
    i = e.agentid
    if resists_infection(i, e, u.d)
        return simEvent[] ## agent not susceptible
    else
        push!(d.I_set, i) ## append to current infected set
        update_counts!(u.d, i)
        duration = rand(d.t_recovery)
        recov_event = becomeRecovered(e.t + duration, i)

        neigh_holiday::Vector{UInt32} = findall(view(d.netw[:holiday],:,i)) 
        ## possible ephemeral/location-based contacts
        
        neigh_loc_holiday::Vector{UInt32} = let loc_idx = loc_holiday(u.d,i); loc_idx > 0 ? findall(view(d.netw[:loc_matrix_holiday],:,loc_idx)) : UInt32[] end

        ## cheat to keep the overall infection rate unchanged
        ## doing this instead of just reading the params because distr fn might not really be a distr fn
        neigh_hh::Vector{UInt32} = findall(view(d.netw[:hh],:,i)) 
        neigh_wp::Vector{UInt32} = findall(view(d.netw[:wp],:,i)) 
        neigh_sch::Vector{UInt32} = findall(view(d.netw[:sch],:,i)) 
        neigh_gq::Vector{UInt32} = findall(view(d.netw[:gq],:,i)) 
        neigh_loc_res::Vector{UInt32} = let loc_idx = loc_res(u.d,i); loc_idx > 0 ? findall(view(d.netw[:loc_matrix_res],:,loc_idx)) : UInt32[] end
        neigh_loc_work::Vector{UInt32} = let loc_idx = loc_work(u.d,i); loc_idx > 0 ? findall(view(d.netw[:loc_matrix_work],:,loc_idx)) : UInt32[] end
        N_contacts_hh = length(d.distr_fns[:hh](neigh_hh, d.distr_params[:hh]))
        N_contacts_wp = length(d.distr_fns[:non_hh](neigh_wp, d.distr_params[:non_hh]))
        N_contacts_sch = length(d.distr_fns[:non_hh](neigh_sch, d.distr_params[:non_hh]))
        N_contacts_gq = length(d.distr_fns[:non_hh](neigh_gq, d.distr_params[:non_hh]))
        N_contacts_loc_res = length(d.distr_fns[:loc_res](neigh_loc_res, d.distr_params[:loc_res]))
        N_contacts_loc_work = length(d.distr_fns[:loc_work](neigh_loc_work, d.distr_params[:loc_work]))

        ## separate samples in case p_inf's are different
        contacts_non_hh::Vector{UInt32} = rsamp(neigh_holiday, N_contacts_wp+N_contacts_sch+N_contacts_gq)
        contacts_hh::Vector{UInt32} = rsamp(neigh_holiday, N_contacts_hh)
        contacts_loc_res::Vector{UInt32} = rsamp(neigh_loc_holiday, N_contacts_loc_res)
        contacts_loc_work::Vector{UInt32} = rsamp(neigh_loc_holiday, N_contacts_loc_work)
        ## unchanged
        #contacts_nonloc::Vector{UInt32} = d.distr_fns[:nonloc](axes(d.netw[:hh],1), d.distr_params[:nonloc]) 

        infected::Vector{UInt32} = unique(vcat(
            randsubseq(contacts_non_hh,d.inf_probs[:p_inf]),
            randsubseq(contacts_hh,d.inf_probs[:p_inf_hh]),
            randsubseq(contacts_loc_res,d.inf_probs[:p_inf_loc]),
            randsubseq(contacts_loc_work,d.inf_probs[:p_inf_loc])#,randsubseq(contacts_nonloc,d.inf_probs[:p_inf_loc])
            ))

        infection_events = [infectionEvent(e.t + rand(0:duration), targ) for targ in infected]
        ##
        ## TODO: log number of secondary infections
        ##

        ## return infection events and recovery event
        return [infection_events; recov_event]
    end
end


function handle_event!(u::SimUnit, e::becomeRecovered)::Vector{simEvent}
    ## change state; currently no future event
    delete!(u.d.I_set, e.agentid) ## deleting from a set is O(1)
    push!(u.d.R_set, e.agentid) ## add to recovered set
    return simEvent[]
end

## periodStuff handles out-of-network infections (people who work or live outside of the synth area)
function handle_event!(u::SimUnit, e::periodicStuff)::Vector{simEvent}

    d = u.d
    ## proportion infected in unit's pop
    ##
    ## TODO: this doesn't actually estimate the proportion of the population infected, does that matter?
    ##
    P_infected = length(d.I_set) / d.n_agents

    ## mean number of infected home / work connections
    
    ##
    ## TODO: add mean location-based contacts to this calculation
    ##

    n_home = P_infected * d.mean_hh_connections
    n_work = P_infected * d.mean_wp_connections

    ## prob of getting infected at home / work for person with unknown home / work connections
    ##  assumes we're checking every (mean contagious duration) days, and p_inf is on that time scale
    ##
    ## TODO: this calculation assumes # of contacts = # of connections; use the distr_fn instead
    ##
    p_not_w = 1.0 - d.inf_probs[:p_inf]
    p_not_h = 1.0 - d.inf_probs[:p_inf_hh]
    p_home = 1.0 - p_not_h^n_home
    p_work = 1.0 - p_not_w^n_work

    ## suspend work infection during holiday time
    #if in(e.t, d.intervals[:holiday])
    #    p_work = 0.0
    #end

    outw_targets = d.outw_assigned
    dumm_targets = d.dummies_assigned

    ## out-workers get infected at non-existent workplace; dummies get infected at nonexistent home 
    infected = [randsubseq(outw_targets, p_work); randsubseq(dumm_targets, p_home)]

    ## time point is random between now and next time this event occurs
    stuff_period = d.periodic_stuff_period
    infection_events = [infectionEvent(e.t + rand(0:stuff_period-1), targ) for targ in infected]

    ## queue for next period
    q_event!(u, periodicStuff(e.t + stuff_period))
    return infection_events
end

## reporting
function handle_event!(u::SimUnit, e::reportingEvent)::Vector{simEvent}
    d = u.d
    ## append current counts to report series
    push!(d.report_series[:active], length(d.I_set))
    push!(d.report_series[:cumI], d.cumI[:agents_assigned])
    push!(d.report_series[:cumI_low_inc_wp], d.cumI[:low_inc_workplace])
    push!(d.report_by_geo, copy(d.cumI_by_geo)) ## copy vector to append current values
    ## queue for next period
    q_event!(u, reportingEvent(e.t + d.report_freq))
    return simEvent[]
end

## tells local processes what to return at the end (returning the whole simunit might take too much memory)
function summarize(u::SimUnit)
    d = u.d
    return Dict(
        :active => d.report_series[:active], 
        :cumI =>  d.report_series[:cumI], 
        :cumI_low_inc_wp => d.report_series[:cumI_low_inc_wp],
        :cumI_by_geo => d.report_by_geo,
        :endR => length(d.R_set),
        :q_len => length(u.q),
        :glob_len => length(u.global_events),
        :n_agents => d.n_agents
    )
end

## event sorting

## periodStuff handler is written to generate events only for local agents,
## so send all events originating from periodStuff to the local queue
##  (dispatching on type of "orig" originating event; subtype dispatch takes precedence over supertype)
function sort_events!(local_sim::SimUnit, orig::periodicStuff, events::Vector{simEvent})
    for e in events
        q_event!(local_sim, e)
    end
    return nothing
end

## infectionEvent generates a single contagious event, which is the same agent and therefore on the same unit
function sort_events!(local_sim::SimUnit, orig::infectionEvent, events::Vector{simEvent})
    for e in events
        q_event!(local_sim, e)
    end
    return nothing
end

## for events generated by other event handlers, place them in local queue, or save for global sync
## depending on criteria defined in this function
## local queue and global cache are inside "local_sim" dict
function sort_events!(local_sim::SimUnit, orig::simEvent, events::Vector{simEvent})

    for e in events
        if e.agentid in local_sim.d.sets[:agents_assigned] ## this check is O(1)
            q_event!(local_sim, e)
        else
            push!(local_sim.global_events, e)
        end
    end

    return nothing
end

end ## @everywhere begin

## only the first worker (global sync overseer) gets the defintions below

## data needed to perform global sync
## using a struct so it will have a strict type
struct globalData <: abstractGlobData
    ## a vector indicating which sim unit owns every agent
    idxs::Vector{UInt8}
end

## globalData constructor, called in gabm.jl spawn_units!()
## should return whatever sort_glob_events() needs
function globalData(inputs)
    g = zeros(UInt8,size(inputs[:netw_hh],2));
    for (k,r) in inputs[:agents_assigned]
        for i in collect(r)
            g[i] = k
        end
    end
    return globalData(g)
end

## direct events from global queue to a local sim unit based on the criteria in this function
## on sync, each sim unit sends a vector of global events; this function is called once for each such vector
## returns a dict specifying which sim unit will get which events
function sort_glob_events(glob_data::globalData, targets::Vector{Int64}, global_events::Vector{simEvent})::Dict{Int64,Vector{simEvent}}
    d = Dict(i => simEvent[] for i in targets) ## empty vectors if no events
    for e in global_events
        unit_id = glob_data.idxs[e.agentid]
        if haskey(d, unit_id)
            push!(d[unit_id], e)
        end
    end
    return d
end

## mean household connections, excluding household-less dummies
calc_mean_hh(netw_mat,exclude) = mean(sum(netw_mat[:,setdiff(axes(netw_mat,2), exclude)], dims=1))
## mean workplace connections = mean non-hh cnxs for people with 1+ non-hh cnxs, excluding workplace-less out-workers
calc_mean_wp(netw_mat,exclude) = mean(filter(x->x>0, sum(netw_mat[:,setdiff(axes(netw_mat,2), exclude)], dims=1)))

## assigning agents to sim units; just random for now
function agent_assignment(unit_ids::Vector{Int64}, nn::Int64)
    n = length(unit_ids)
    ## not really necessary to randomize but doesn't hurt
    splits = ranges(lrRound(fill(nn/n, n)))
    idxs = UInt32.(shuffle(1:nn))
    return Dict(unit_ids .=> [idxs[i] for i in splits])
end

## data needed to initialize sim units;  called in gabm.jl spawn_units!()
## all external inputs go into this dict, which the global process then uses to initialize local sim units
## (this may seem a little convoluted, but it allows the main process to read in whatever data is needed,
##    copy pieces to the local workers, and then delete this object to free up memory)
function modelInputs(unit_ids::Vector{Int64}; kwargs...)
    ## note, calling kwargs like a fn is enabled by a definition in utils.jl
    netw_hh = SparseMatrixCSC{Bool, UInt32}(kwargs(:netw_hh, ()->Symmetric(dser_path("jlse/adj_mat_hh.jlse"))))
    netw_wp = SparseMatrixCSC{Bool, UInt32}(kwargs(:netw_wp, ()->Symmetric(dser_path("jlse/adj_mat_wp.jlse"))))
    ## dummies appear in work networks, but live outside the synth pop (no household or demographic info generated)
    ## local sim unit is responsible for determining if they got infected at home
    dummies = Set{UInt32}(kwargs(:dummies, ()->keys(dser_path("jlse/adj_dummy_keys.jlse"))))
    ## outside workers have households, but no workplace network
    ## local sim unit is responsible for determining if they got infected at work
    outw = Set{UInt32}(kwargs(:out_workers, ()->keys(dser_path("jlse/adj_out_workers.jlse"))))

    println("returning model inputs")
    ## this constructor syntax is enabled by @kwdef macro above
    return Dict(
        :netw_hh => netw_hh,
        :netw_wp => netw_wp,
        :netw_sch => SparseMatrixCSC{Bool, UInt32}(kwargs(:netw_sch, ()->Symmetric(dser_path("jlse/adj_mat_sch.jlse")))),
        :netw_gq => SparseMatrixCSC{Bool, UInt32}(kwargs(:netw_gq, ()->Symmetric(dser_path("jlse/adj_mat_gq.jlse")))),
        :loc_matrix_res => SparseMatrixCSC{Bool,UInt32}(kwargs(:loc_matrix_res, ()->dser_path("jlse/res_loc_contact_mat.jlse"))),
        :loc_matrix_work => SparseMatrixCSC{Bool,UInt32}(kwargs(:loc_matrix_work, ()->dser_path("jlse/work_loc_contact_mat.jlse"))),
        :loc_lookup_res => Dict{UInt32,UInt32}(kwargs(:loc_lookup_res, ()->dser_path("jlse/res_loc_lookup.jlse"))),
        :loc_lookup_work => Dict{UInt32,UInt32}(kwargs(:loc_lookup_work, ()->dser_path("jlse/work_loc_lookup.jlse"))),
        :geo_lookup => Dict{UInt32,UInt16}(kwargs(:geo_lookup, ()->dser_path("precalc/cbg_idx_lookup.jlse"))),
        :t_inc => UInt32(get(kwargs, :t_inc, 5)),
        :t_recovery => get(kwargs, :t_recovery, 8:12),
        :init_inf => get(kwargs, :init_inf, [first(unit_ids) => 10]),
        :agents_assigned => kwargs(:agents_assigned, ()->agent_assignment(unit_ids, size(netw_hh,2))),
        :dummies => dummies,
        :out_workers => outw,
        :commuters_to_low_inc_wp => Set{UInt32}(kwargs(:low_inc_workplace, ()->dser_path("precalc/low_inc_workplace.jlse"))),
        :k12_worker => Set{UInt32}(kwargs(:k12_worker, ()->dser_path("precalc/k12_worker.jlse"))),
        :mean_hh_connections => Float64(kwargs(:mean_hh_connections, ()->calc_mean_hh(netw_hh,dummies))),
        :mean_wp_connections => Float64(kwargs(:mean_wp_connections, ()->calc_mean_wp(netw_wp,outw))),
        :report_freq => UInt32(get(kwargs, :report_freq, 5)),
        ## defaults for within-household infection
        :p_inf_hh => get(kwargs, :p_inf_hh, 0.15), 
        :distr_fn_hh => get(kwargs, :distr_fn_hh, :const),
        :distr_params_hh => get(kwargs, :distr_params_hh, (16,)),
        ## defaults for work,school,GQ infection
        :p_inf => get(kwargs, :p_inf, 0.15),
        :distr_fn_non_hh => get(kwargs, :distr_fn_non_hh, :const),
        :distr_params_non_hh => get(kwargs, :distr_params_non_hh, (8,)),
        ## defaults for ephemeral/location-based infection
        :p_inf_loc => get(kwargs, :p_inf_loc, 0.15),
        :distr_fn_loc_res => get(kwargs, :distr_fn_loc_res, :zero),
        :distr_params_loc_res => get(kwargs, :distr_params_loc_res, ()),
        :distr_fn_loc_work => get(kwargs, :distr_fn_loc_work, :zero),
        :distr_params_loc_work => get(kwargs, :distr_params_loc_work, ()),
        :distr_fn_nonloc => get(kwargs, :distr_fn_nonloc, :zero),
        :distr_params_nonloc => get(kwargs, :distr_params_nonloc, ()),
        ## special time intervals; missing by default
        :nonessential_wp_closed => get(kwargs, :nonessential_wp_closed, missing),
        :sch_closed => get(kwargs, :sch_closed, missing),
        ## misc flags
        :flags => Set(get(kwargs, :flags, Symbol[]))
        )
end

## queue initial infections when init_inf is [unit id => # infections]
function q_init_inf!(u::SimUnit, init_inf::Vector{Pair{I,I}}) where I<:Integer
    id = u.id
    d_init = Dict(init_inf)
    if haskey(d_init, id)
        for i in rand(collect(u.d.sets[:agents_assigned]), d_init[id])
            #q_event!(u, infectionEvent(1,i))
            q_event!(u, becomeContagious(0,i))
        end
    end
end

## queue initial infections when init_inf is a list of agent ids
function q_init_inf!(u::SimUnit, init_inf::Vector{I}) where I<:Integer
    for i in init_inf
        if i in u.d.sets[:agents_assigned]
            #q_event!(u, infectionEvent(1,i))
            q_event!(u, becomeContagious(0,i))
        end
    end
end

## a subset of a sparse matrix, but keep size and indices unchanged
## (everything else becomes zeros, which use no memory)
function sp_col_subset(orig::SparseMatrixCSC{T,U}, col_idxs) where {T<:Any,U<:Integer}
    r = spzeros(T,U,size(orig))
    r[:,col_idxs] = copy(orig[:,col_idxs])
    return r
end

## this fn is called once for each sim unit, from a fn in gabm.jl
## adds domain-specific data to a sim unit
## adds initial event(s) to queue (at least one unit must have an initial event)
function init_sim_unit!(u::SimUnit, inputs::Dict{Symbol,Any})

    id = u.id
    print("init sim unit ", id, "; ")
    u.t_inc = inputs[:t_inc] ## defined for SimUnit not SimData, because global process needs it

    d = SimData(); u.d = d ## make the struct and give sim unit the reference
    agents_assigned = inputs[:agents_assigned][id]

    d.n_agents = length(agents_assigned)
    d.I_set = Set{UInt32}()
    d.R_set = Set{UInt32}()
    d.t_recovery = inputs[:t_recovery]
    ## infection probability is defined in terms of infectiousness duration, so 
    ##  update out-of-network infections on the same timescale so probabilities are correct
    d.periodic_stuff_period = round(UInt32, mean(inputs[:t_recovery]))
    d.report_freq = inputs[:report_freq]

    ## having many network matrices wastes memory because indices are repeated; could combine into an 8-bit binary
    d.netw = Dict(
        :hh => sp_col_subset(inputs[:netw_hh], agents_assigned),
        :wp => sp_col_subset(inputs[:netw_wp], agents_assigned),
        :sch => sp_col_subset(inputs[:netw_sch], agents_assigned),
        :gq => sp_col_subset(inputs[:netw_gq], agents_assigned),
        ### copy() will happen automatically anyway?
        :loc_matrix_res => copy(inputs[:loc_matrix_res]), ## local ephemeral contacts by location; every sim unit needs the full matrix
        :loc_matrix_work => copy(inputs[:loc_matrix_work]),
        ## non-local contact network is everyone in the population, so just a single column of "trues" 
        ## TODO: this wastes memory; but we're not currently using it
        ## :nonlocal => SparseMatrixCSC{Bool, UInt32}(trues(size(inputs[:netw_hh], 1), 1))
    )

    ## agent assignment and other membership groups
    ##   if there are many of these, could change representation to bitmatrix or something 
    d.sets = Dict(
        :agents_assigned => Set{UInt32}(agents_assigned), ## must be something wih O(1) lookup
        :low_inc_workplace => intersect(inputs[:commuters_to_low_inc_wp], agents_assigned),
        :k12_worker => intersect(inputs[:k12_worker] , agents_assigned)
    )

    d.inf_probs = Dict(
        :p_inf => inputs[:p_inf],
        :p_inf_hh => inputs[:p_inf_hh],
        :p_inf_loc => inputs[:p_inf_loc]
    )

    ## note, loc_lookup_res is missing inst gq residents (and dummies); loc_lookup_work has only commuters
    d.loc_lookup = Dict(
        :res => dsubset(inputs[:loc_lookup_res], agents_assigned), ## don't need everyone's lookup, save some memory
        :work => dsubset(inputs[:loc_lookup_work], agents_assigned)
    )

    ## which function to use to determine contact events
    opts = Dict(:all=>distr_all, :const=>distr_const, :zero=>distr_zero, :pois=>distr_pois, :exp=>distr_exp)

    d.distr_fns = Dict(
        :hh => get(opts, inputs[:distr_fn_hh], distr_const),
        :non_hh => get(opts, inputs[:distr_fn_non_hh], distr_const),
        :loc_res => get(opts, inputs[:distr_fn_loc_res], distr_zero),
        :loc_work => get(opts, inputs[:distr_fn_loc_work], distr_zero),
        :nonloc => distr_zero ## currently off; get(opts, inputs[:distr_fn_nonloc], distr_zero)
    )

    d.distr_params = Dict(
        :hh => inputs[:distr_params_hh],
        :non_hh => inputs[:distr_params_non_hh],
        :loc_res => inputs[:distr_params_loc_res],
        :loc_work => inputs[:distr_params_loc_work],
        :nonloc => inputs[:distr_params_nonloc]
    )

    d.dummies_assigned = collect(UInt32, intersect(inputs[:dummies], agents_assigned)) ## collect() because randsubseq can't handle Set
    d.outw_assigned = collect(UInt32, intersect(inputs[:out_workers], agents_assigned)) ## collect() because randsubseq can't handle Set
    ## note, these are currently global means (not per sim unit)
    d.mean_hh_connections = inputs[:mean_hh_connections]
    d.mean_wp_connections = inputs[:mean_wp_connections]

    ## groups for which data should be reported
    d.report_groups = collect(keys(d.sets)) ## all groups by default
    ## report cumulative infections for those groups
    d.cumI = Dict(d.report_groups .=> 0)
    ## initialize time series, modified in reporting event handler
    d.report_series = Dict(k => Int[] for k in [:active, :cumI, :cumI_low_inc_wp])

    ## for tracking # infections by home location
    d.geo_assigned = dsubset(inputs[:geo_lookup], agents_assigned)
    d.cumI_by_geo = zeros(Int, maximum(values(inputs[:geo_lookup])))
    d.report_by_geo = Vector{Int}[]

    ## include keys only for specified intervals
    d.intervals = Dict{Symbol, UnitRange{Int64}}(k=>inputs[k] for k in 
        [:nonessential_wp_closed, :sch_closed] if !ismissing(inputs[k]))
    
    d.flags = inputs[:flags]

    println("using contact functions ", d.distr_fns[:hh], " ", d.distr_fns[:non_hh], " ", d.distr_fns[:loc_res], " ", d.distr_fns[:loc_work], " ", d.distr_fns[:nonloc])
    println("using contact params ", d.distr_params[:hh], " ", d.distr_params[:non_hh], " ", d.distr_params[:loc_res], " ", d.distr_params[:loc_work], " ", d.distr_params[:nonloc], " ", d.flags)

    ## queue initial events
    q_init_inf!(u, inputs[:init_inf])
    q_event!(u, periodicStuff(d.periodic_stuff_period))
    q_event!(u, reportingEvent(1)) ## offset from sync

    return nothing
end


## pre-calculate stuff used in intialization

function precalc_geo()
    p_idxs = let k = dser_path("jlse/adj_mat_keys.jlse"); Dict(k .=> UInt32.(eachindex(k))); end
    d = Dict(i=>k[3] for (k,i) in p_idxs)
    ser_path("precalc/cbg_idx_lookup.jlse",d)
    return nothing
end

## can store people in low-inc worplaces, schools, etc., as a set of network indices
##   if using several of these sets, could save memory (with only a small performace cost) by combining
##   them with agents_assigned into a single dict
function precalc_sets()
    p_idxs = let k = dser_path("jlse/adj_mat_keys.jlse"); Dict(k .=> UInt32.(eachindex(k))); end
    wp_low_inc = let w = dser_path("jlse/company_workers.jlse"); filterk(x->x[2]==UInt8(1), w); end
    people = dser_path("jlse/people.jlse")
    work_dummies = dser_path("jlse/work_dummies.jlse")
    sch_workers = dser_path("jlse/sch_workers.jlse")

    println("calc stats")
    ## find person keys in each category
    commute_to_low_inc_wp = Set(p_idxs[k[1:3]] for k in reduce(vcat, collect(values(wp_low_inc))))

    low_inc_and_commute = union(Set(p_idxs[k[1:3]] for k in keys(filterv(x->x.com_LODES_low, people))),
        Set(p_idxs[k[1:3]] for k in filter(x->x[4]==1, work_dummies)))

    p_in_school = Set(p_idxs[k[1:3]] for k in keys(filterv(p->(!ismissing(p.sch_grade) && !in(p.sch_grade, ["c","g"])), people)))
    work_in_school = Set(p_idxs[k[1:3]] for k in reduce(vcat, collect(values(sch_workers))))
    school_student_or_worker = union(p_in_school,work_in_school)
    age65o = Set(p_idxs[k[1:3]] for k in keys(filterv(x->x.age>=65, people)))

    ser_path("precalc/low_inc_workplace.jlse", commute_to_low_inc_wp)
    ser_path("precalc/low_inc_commuter.jlse", low_inc_and_commute)
    ser_path("precalc/k12_student.jlse", p_in_school)
    ser_path("precalc/k12_worker.jlse", work_in_school)
    ser_path("precalc/k12_student_or_worker.jlse", school_student_or_worker)
    ser_path("precalc/age_65o.jlse", age65o)    
    return nothing
end

function precalc_stats()
    p_idxs = let k = dser_path("jlse/adj_mat_keys.jlse"); Dict(k .=> eachindex(k)); end
    wp_low_inc = let w = dser_path("jlse/company_workers.jlse"); filterk(x->x[2]==UInt8(1), w); end
    people = dser_path("jlse/people.jlse")
    work_dummies = dser_path("jlse/work_dummies.jlse")
    sch_workers = dser_path("jlse/sch_workers.jlse")

    println("calc stats")
    ## find person keys in each category
    commute_to_low_inc_wp = Set(k[1:3] for k in reduce(vcat, collect(values(wp_low_inc))))

    low_inc_and_commute = union(Set(keys(filterv(x->x.com_LODES_low, people))),
        Set(k[1:3] for k in filter(x->x[4]==1, work_dummies)))

    p_in_school = Set(keys(filterv(p->(!ismissing(p.sch_grade) && !in(p.sch_grade, ["c","g"])), people)))
    work_in_school = Set(k[1:3] for k in reduce(vcat, collect(values(sch_workers))))
    school_student_or_worker = union(p_in_school,work_in_school)
    age65o = Set(keys(filterv(x->x.age>=65, people)))

    ## boolean vector of category memberships for each person, stored as an 8-bit integer (to conserve memory)
    ## the keys of this dict are each person's matrix index
    statdict = Dict(UInt32(i) => UInt8(evalpoly(2,[in(k,s) for s in 
        [keys(p_idxs), commute_to_low_inc_wp, low_inc_and_commute, p_in_school, school_student_or_worker, age65o]]))
        for (k,i) in p_idxs)

    ## store category names in the same order
    statkeys = let x = [:exists, :low_inc_workplace, :low_inc_commuter, :k12_student, :k12_student_or_worker, :age_65o];
        Dict(x .=> eachindex(x)) end

    ser_path("precalc/person_stats.jlse", statdict)
    ser_path("precalc/stat_keys.jlse",statkeys)
    return nothing
end

## use bit arithmetic to check membership in category name k
function check_stat(person_stats::Dict{UInt32,UInt8}, stat_keys::Dict{Symbol,Int}, person_idx::UInt32, k::Symbol)
    return person_stats[person_idx] & (UInt8(1) << stat_keys[k] >> 1) != 0
end




