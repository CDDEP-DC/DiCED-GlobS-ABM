using LinearAlgebra
include("utils.jl")
include("fileutils.jl")

## remote workers need these definitions
##
## NOTE : remote workers currently don't have utils.jl (to save memory)
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

##
## event handlers, dispatch on event type
## worker loop expects these functions to return future events as a Vector{simEvent}
##

## most efficient way to do logging/reporting:
## add code to existing events (e.g., increment case count in E->I event)
##  and queue a "data summary" event similar to sync event

function handle_event!(local_sim::SimUnit, e::infectionEvent)::Vector{simEvent}
    ## target may be in another sim unit, so no state change yet
    ## just return an event for the target to become contagious at the correct time
    ## (note if an agent becomes infected twice, the earlier one will take effect & the later will do nothing)
    return [becomeContagious(e.t + local_sim[:t_inc], e.agentid)]
end

## no longer doing it this way
## Bernouilli sampling, everyone in "neigh" has probability p of becoming infected
#function infect_Bern(neigh::Vector{UInt32}, p::Float64)::Vector{UInt32}
#    return randsubseq(neigh, p) ## fn provided by Random, promises efficient Bernouilli sampling
#end
## people encounter some contacts often, others not at all
## overall transmission prob same as bern, but small # of transmissions more likely
#function infect_clumped(neigh::Vector{UInt32}, p::Float64)::Vector{UInt32}
#    return unique(randsubseq(rand(neigh,length(neigh)), p)) ## rand() samples with replacement
#end

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

##
## the becomeContagious handler
##   if agent is susceptible, change state, and generate all future infection events that result
##
function handle_event!(u::SimUnit, e::becomeContagious)::Vector{simEvent}
    i = e.agentid
    if in(i, u[:I_set]) || in(i, u[:R_set]) ## testing set membership is O(1)
        return simEvent[] ## agent not susceptible
    else
        u[:cum_I] += 1 ## update cumulative infection count
        push!(u[:I_set], i) ## append to current infected set
        duration = rand(u[:t_recovery])
        recov_event = becomeRecovered(e.t + duration, i)

        neigh_non_hh::Vector{UInt32} = findall(u[:netw_non_hh][:,i]) ## note, this is not really a "find", it's just a lookup in a sparse array
        neigh_hh::Vector{UInt32} = findall(u[:netw_hh][:,i]) 

        ## possible ephemeral/location-based contacts
        if haskey(u[:loc_lookup_res], i)
            res_loc::UInt32 = u[:loc_lookup_res][i]
            neigh_loc_res::Vector{UInt32} = findall(u[:loc_matrix_res][:,res_loc]) ## columns of this sparse matrix are locations
        else
            neigh_loc_res = UInt32[]
        end
        if haskey(u[:loc_lookup_work], i)
            work_loc::UInt32 = u[:loc_lookup_work][i]
            neigh_loc_work::Vector{UInt32} = findall(u[:loc_matrix_work][:,work_loc])
        else
            neigh_loc_work = UInt32[]
        end

        contacts_non_hh::Vector{UInt32} = u[:distr_fn_non_hh](neigh_non_hh, u[:distr_params_non_hh])
        contacts_hh::Vector{UInt32} = u[:distr_fn_hh](neigh_hh, u[:distr_params_hh])

        ## if keeping overall contact rate constant, don't allow mean loc contacts to be greater than # neighbors
        if u[:h_test] && (get(u[:distr_params_loc_res],1,0) > length(neigh_hh))
            distr_params_loc_res = (length(neigh_hh), u[:distr_params_loc_res][2:end]...)
        else
            distr_params_loc_res = u[:distr_params_loc_res]
        end

        if u[:w_test] && (get(u[:distr_params_loc_work],1,0) > length(neigh_non_hh))
            distr_params_loc_work = (length(neigh_non_hh), u[:distr_params_loc_work][2:end]...)
        else
            distr_params_loc_work = u[:distr_params_loc_work]
        end

        if u[:h_test] && (get(u[:distr_params_nonloc],1,0) > length(neigh_hh))
            distr_params_nonloc = (length(neigh_hh), u[:distr_params_nonloc][2:end]...)
        elseif u[:w_test] && (get(u[:distr_params_nonloc],1,0) > length(neigh_non_hh))
            distr_params_nonloc = (length(neigh_non_hh), u[:distr_params_nonloc][2:end]...)
        else
            distr_params_nonloc = u[:distr_params_nonloc]
        end

        contacts_loc_res::Vector{UInt32} = u[:distr_fn_loc_res](neigh_loc_res, distr_params_loc_res)
        contacts_loc_work::Vector{UInt32} = u[:distr_fn_loc_work](neigh_loc_work, distr_params_loc_work)
        ## anyone in the population is a potential nonlocal contact
        contacts_nonloc::Vector{UInt32} = u[:distr_fn_nonloc](axes(u[:netw_non_hh],1), distr_params_nonloc) 
        ## randsubseq promises efficient Bernouilli sampling
        ##   unique() because infecting someone twice doesn't have any additional effect
        infected::Vector{UInt32} = unique(vcat(
            randsubseq(contacts_non_hh,u[:p_inf]),
            randsubseq(contacts_hh,u[:p_inf_hh]),
            randsubseq(contacts_loc_res,u[:p_inf_loc]),
            randsubseq(contacts_loc_work,u[:p_inf_loc]),
            randsubseq(contacts_nonloc,u[:p_inf_loc])))

        infection_events = [infectionEvent(e.t + rand(0:duration), targ) for targ in infected]
        ##
        ## TODO: log number of secondary infections
        ##

        ## return infection events and recovery event
        return [infection_events; recov_event]
    end

end

function handle_event!(local_sim::SimUnit, e::becomeRecovered)::Vector{simEvent}
    ## change state; currently no future event
    delete!(local_sim[:I_set], e.agentid) ## deleting from a set is O(1)
    push!(local_sim[:R_set], e.agentid) ## add to recovered set
    return simEvent[]
end

## periodStuff handles out-of-network infections (people who work or live outside of the synth area)
function handle_event!(u::SimUnit, e::periodicStuff)::Vector{simEvent}

    ## proportion infected in unit's pop
    ##
    ## TODO: this doesn't actually estimate the proportion of the population infected, does that matter?
    ##
    P_infected = length(u[:I_set]) / u[:n_agents]

    ## mean number of infected home / work connections
    
    ##
    ## TODO: add mean location-based contacts to this calculation
    ##

    n_home = P_infected * u[:mean_hh_connections]
    n_work = P_infected * u[:mean_wp_connections]

    ## prob of getting infected at home / work for person with unknown home / work connections
    ##  assumes we're checking every (mean contagious duration) days, and p_inf is on that time scale
    p_not = 1.0 - u[:p_inf]
    p_home = 1.0 - p_not^n_home
    p_work = 1.0 - p_not^n_work

    ## out-workers get infected at non-existent workplace; dummies get infected at nonexistent home 
    infected = [randsubseq(u[:outw_assigned], p_work); 
                randsubseq(u[:dummies_assigned], p_home)]
 
    ## time point is random between now and next time this event occurs
    infection_events = [infectionEvent(e.t + rand(0:u[:periodic_stuff_period]-1), targ) for targ in infected]

    ## queue for next period
    q_event!(u, periodicStuff(e.t + u[:periodic_stuff_period]))
    return infection_events
end

## reporting
function handle_event!(u::SimUnit, e::reportingEvent)::Vector{simEvent}

    ## append number infected to log
    push!(u[:log], (e.t, length(u[:I_set])))

    ## queue for next period
    q_event!(u, reportingEvent(e.t + u[:report_freq]))
    return simEvent[]
end

## periodStuff handler is written to generate events only for local agents,
## so send all events originating from periodStuff to the local queue
##  (dispatching on type of "orig" originating event; subtype dispatch takes precedence over supertype)
function sort_events!(local_sim::SimUnit, orig::periodicStuff, events::Vector{simEvent})
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
        if e.agentid in local_sim[:agents_assigned] ## this check is O(1)
            q_event!(local_sim, e)
        else
            push!(local_sim[:global_events], e)
        end
    end

    return nothing
end

## tells local processes what to return at the end (returning the whole simunit might take too much memory)
function summarize(u::SimUnit)
    return Dict(
        :cum_I => u[:cum_I],
        :curr_I => length(u[:I_set]),
        :curr_R => length(u[:R_set]),
        :q_len => length(u[:q]),
        :glob_len => length(u[:global_events]),
        :n_agents => u[:n_agents],
        :log => deepcopy(u[:log])
    )
end

end ## @everywhere begin
## only the first worker (global sync overseer) gets the defintions below

## data needed to perform global sync
## using a struct so it will have a strict type
struct globalData <: abstractGlobData
    ## a vector indicating which sim unit owns every agent
    idxs::Vector{UInt8}
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

## data needed to initialize sim units
## all external inputs go into this struct, which the global process then uses to initialize local sim units
##   parametric types T U V etc can be anything the functions below are written to handle
##   (this didn't really need to be a struct, but it makes the syntax below cleaner)
@kwdef struct modelInputs{U,V,W,Ta,Tb,Tc,Td,Te} <: abstractModelInputs
    #netw::SparseMatrixCSC{Bool,UInt32}
    netw_hh::SparseMatrixCSC{Bool,UInt32}
    netw_non_hh::SparseMatrixCSC{Bool,UInt32}
    loc_matrix_res::SparseMatrixCSC{Bool,UInt32}
    loc_matrix_work::SparseMatrixCSC{Bool,UInt32}
    loc_lookup_res::Dict{UInt32,UInt32}
    loc_lookup_work::Dict{UInt32,UInt32}
    agents_assigned::U
    t_inc::Int64
    p_inf::Float64
    p_inf_hh::Float64
    p_inf_loc::Float64
    t_recovery::V
    init_inf::Dict{Int64,Int64}
    dummies_assigned::W
    outw_assigned::W
    mean_hh_connections::Float64
    mean_wp_connections::Float64
    report_freq::Int64
    #infection_function::Symbol
    distr_fn_hh::Symbol
    distr_params_hh::Ta
    distr_fn_non_hh::Symbol
    distr_params_non_hh::Tb
    distr_fn_loc_res::Symbol
    distr_params_loc_res::Tc
    distr_fn_loc_work::Symbol
    distr_params_loc_work::Td
    distr_fn_nonloc::Symbol
    distr_params_nonloc::Te
    flags::Set{Symbol}
end

## constructor, called in gabm.jl spawn_units!()
function modelInputs(unit_ids::Vector{Int64}; kwargs...)

    ## these funcs only evaluated when needed
    read_dummies() = Set{UInt32}(keys(dser_path("jlse/adj_dummy_keys.jlse")))
    read_outw() = Set{UInt32}(keys(dser_path("jlse/adj_out_workers.jlse")))
    read_netw_hh() = (println("read hh net"); SparseMatrixCSC{Bool,UInt32}(Symmetric(dser_path("jlse/hh_adj_mat.jlse"))) )
    read_netw_non_hh() = (println("read non-hh net"); SparseMatrixCSC{Bool,UInt32}(Symmetric(dser_path("jlse/adj_mat.jlse"))) )
    read_loc_matrix_res() = (println("read res loc mat"); SparseMatrixCSC{Bool,UInt32}(dser_path("jlse/res_loc_contact_mat.jlse")) )
    read_loc_matrix_work() = (println("read work loc mat"); SparseMatrixCSC{Bool,UInt32}(dser_path("jlse/work_loc_contact_mat.jlse")) )
    read_loc_lookup_res() = Dict{UInt32,UInt32}(dser_path("jlse/res_loc_lookup.jlse"))
    read_loc_lookup_work() = Dict{UInt32,UInt32}(dser_path("jlse/work_loc_lookup.jlse"))

    ## mean household connections, excluding household-less dummies
    calc_mean_hh(netw_hh,dummies) = (println("calc mean hh"); mean(sum(netw_hh[:,setdiff(axes(netw_hh,2), dummies)], dims=1)))
    ## mean workplace connections = mean non-hh cnxs for people with 1+ non-hh cnxs, excluding workplace-less out-workers
    calc_mean_wp(netw_non_hh,outw) = (println("calc mean wp"); mean(filter(x->x>0, sum(netw_non_hh[:,setdiff(axes(netw_non_hh,2), outw)], dims=1))))
    ## for now, just join household and non-hh networks; keyword arg takes precedence if present
    ## calc_full_net(netw_hh,netw_non_hh)::SparseMatrixCSC{Bool, UInt32} = (println("calc full net"); netw_hh .| netw_non_hh)

    ## awkward syntax but avoids costly evaluation when kwargs has the key
    netw_hh::SparseMatrixCSC{Bool, UInt32} = get(read_netw_hh, kwargs, :netw_hh)
    netw_non_hh::SparseMatrixCSC{Bool, UInt32} = get(read_netw_non_hh, kwargs, :netw_non_hh)
    loc_matrix_res::SparseMatrixCSC{Bool, UInt32} = get(read_loc_matrix_res, kwargs, :loc_matrix_res)
    loc_matrix_work::SparseMatrixCSC{Bool, UInt32} = get(read_loc_matrix_work, kwargs, :loc_matrix_work)
    loc_lookup_res::Dict{UInt32,UInt32} = get(read_loc_lookup_res, kwargs, :loc_lookup_res)
    loc_lookup_work::Dict{UInt32,UInt32} = get(read_loc_lookup_work, kwargs, :loc_lookup_work)

    ## dummies appear in work networks, but live outside the synth pop (no household or demographic info generated)
    ## local sim unit is responsible for determining if they got infected at home
    dummies::Set{UInt32} = get(read_dummies, kwargs, :dummies)
    ## outside workers have households, but no workplace network
    ## local sim unit is responsible for determining if they got infected at work
    outw::Set{UInt32} = get(read_outw, kwargs, :out_workers)

    mu_hh_cnx::Float64 = get(()->calc_mean_hh(netw_hh,dummies), kwargs, :mean_hh_connections)
    mu_work_cnx::Float64 = get(()->calc_mean_wp(netw_non_hh,outw), kwargs, :mean_wp_connections)

    ## assigning agents to sim units
    if haskey(kwargs, :agents_assigned)
        println("assignment ", [(k,length(v)) for (k,v) in kwargs[:agents_assigned]])
        assign_idxs = kwargs[:agents_assigned]
    else
        println("random assignment")
        ## TODO: optimize
        n = length(unit_ids)
        nn = size(netw_hh,2)
        #assign_idxs = Dict(unit_ids .=> ranges(lrRound(fill(nn/n, n))))
        ## not really necessary to randomize but doesn't hurt
        splits = ranges(lrRound(fill(nn/n, n)))
        idxs = UInt32.(shuffle(1:nn))
        assign_idxs = Dict(unit_ids .=> [idxs[i] for i in splits])
    end

    assign_dummies = Dict(i => collect(intersect(dummies, assign_idxs[i])) for i in unit_ids)
    assign_outw = Dict(i => collect(intersect(outw, assign_idxs[i])) for i in unit_ids)

    println("returning model inputs")
    ## this constructor syntax is enabled by @kwdef macro above
    return modelInputs(
        netw_hh = netw_hh,
        netw_non_hh = netw_non_hh,
        loc_matrix_res = loc_matrix_res,
        loc_matrix_work = loc_matrix_work,
        loc_lookup_res = loc_lookup_res,
        loc_lookup_work = loc_lookup_work,    
        #netw = get(()->calc_full_net(netw_hh,netw_non_hh), kwargs, :full_net),
        t_inc = get(kwargs, :t_inc, 5),
        t_recovery = get(kwargs, :t_recovery, 8:12),
        init_inf = Dict(get(kwargs, :init_inf, [first(unit_ids) => 10])),
        agents_assigned = assign_idxs,
        dummies_assigned = assign_dummies,
        outw_assigned = assign_outw,
        mean_hh_connections = mu_hh_cnx,
        mean_wp_connections = mu_work_cnx,
        report_freq = get(kwargs, :report_freq, 5),

        ## defaults for within-household infection
        p_inf_hh = get(kwargs, :p_inf_hh, 0.15), 
        distr_fn_hh = get(kwargs, :distr_fn_hh, :const),
        distr_params_hh = get(kwargs, :distr_params_hh, (16,)),

        ## defaults for work,school,GQ infection
        p_inf = get(kwargs, :p_inf, 0.15),
        distr_fn_non_hh = get(kwargs, :distr_fn_non_hh, :const),
        distr_params_non_hh = get(kwargs, :distr_params_non_hh, (8,)),

        ## defaults for ephemeral/location-based infection
        p_inf_loc = get(kwargs, :p_inf_loc, 0.15),
        distr_fn_loc_res = get(kwargs, :distr_fn_loc_res, :zero),
        distr_params_loc_res = get(kwargs, :distr_params_loc_res, ()),
        distr_fn_loc_work = get(kwargs, :distr_fn_loc_work, :zero),
        distr_params_loc_work = get(kwargs, :distr_params_loc_work, ()),
        distr_fn_nonloc = get(kwargs, :distr_fn_nonloc, :zero),
        distr_params_nonloc = get(kwargs, :distr_params_nonloc, ()),

        ## misc flags
        flags = Set(get(kwargs, :flags, Symbol[]))
        )
end

## globalData constructor, called in gabm.jl spawn_units!()
## should return whatever sort_glob_events() needs
function globalData(inputs::modelInputs)
    g = zeros(UInt8,size(inputs.netw_hh,2));
    for (k,r) in inputs.agents_assigned
        for i in collect(r)
            g[i] = k
        end
    end
    return globalData(g)
end

## this fn is called from a fn in gabm.jl
## adds domain-specific data to a sim unit
## adds initial event(s) to queue (at least one unit must have an initial event)
## (note, SimUnit type is currently just a Dict with symbol keys, can add any key needed by fns above)
function init_sim_unit!(u::SimUnit, id::Int64, inputs::modelInputs)

    print("init sim unit ", id, "; ")
    ## add state variables and parameters
    #u[:netw] = spzeros(Bool,UInt32,size(inputs.netw)) ## zeros use no memory in a sparse array
    #u[:netw][:,inputs.agents_assigned[id]] = copy(inputs.netw[:,inputs.agents_assigned[id]]) ## copy just the columns assigned to this unit

    u[:netw_hh] = spzeros(Bool,UInt32,size(inputs.netw_hh)) ## zeros use no memory in a sparse array
    u[:netw_non_hh] = spzeros(Bool,UInt32,size(inputs.netw_non_hh)) ## zeros use no memory in a sparse array
    u[:netw_hh][:,inputs.agents_assigned[id]] = copy(inputs.netw_hh[:,inputs.agents_assigned[id]]) ## copy just the columns assigned to this unit
    u[:netw_non_hh][:,inputs.agents_assigned[id]] = copy(inputs.netw_non_hh[:,inputs.agents_assigned[id]]) ## copy just the columns assigned to this unit
    
    ## every sim unit needs these
    u[:loc_matrix_res] = inputs.loc_matrix_res
    u[:loc_matrix_work] = inputs.loc_matrix_work
    u[:loc_lookup_res] = inputs.loc_lookup_res
    u[:loc_lookup_work] = inputs.loc_lookup_work

    u[:agents_assigned] = Set{UInt32}(inputs.agents_assigned[id]) ## should be something wih O(1) lookup, like Set or UnitRange
    u[:n_agents] = length(u[:agents_assigned])
    u[:dummies_assigned] = Vector{UInt32}(inputs.dummies_assigned[id]) ## vector because randsubseq can't handle Set
    u[:outw_assigned] = Vector{UInt32}(inputs.outw_assigned[id]) ## vector because randsubseq can't handle Set
    u[:cum_I] = 0
    u[:I_set] = Set{UInt32}()
    u[:R_set] = Set{UInt32}()
    u[:t_inc] = inputs.t_inc
    u[:t_recovery] = inputs.t_recovery
    ## note, these are currently global means (not per sim unit)
    u[:mean_hh_connections] = inputs.mean_hh_connections
    u[:mean_wp_connections] = inputs.mean_wp_connections
    ## infection probability is defined in terms of infectiousness duration, so 
    ##  update out-of-network infections on the same timescale so probabilities are correct
    u[:periodic_stuff_period] = round(Int, mean(inputs.t_recovery))
    u[:report_freq] = inputs.report_freq

    u[:p_inf] = inputs.p_inf
    u[:p_inf_hh] = inputs.p_inf_hh
    u[:p_inf_loc] = inputs.p_inf_loc

    ## which function to use for infections
    ## no longer doing it this way
    #opts = Dict(:bern=>infect_Bern,:clumped=>infect_clumped)
    #u[:infect_fn] = get(opts, inputs.infection_function, infect_Bern)
    #println("using infection function ", u[:infect_fn])

    ## which function to use to determine contact events
    opts = Dict(:all=>distr_all, :const=>distr_const, :zero=>distr_zero, :pois=>distr_pois, :exp=>distr_exp)
    u[:distr_fn_hh] = get(opts, inputs.distr_fn_hh, distr_const)
    u[:distr_params_hh] = inputs.distr_params_hh
    u[:distr_fn_non_hh] = get(opts, inputs.distr_fn_non_hh, distr_const)
    u[:distr_params_non_hh] = inputs.distr_params_non_hh
    u[:distr_fn_loc_res] = get(opts, inputs.distr_fn_loc_res, distr_zero)
    u[:distr_params_loc_res] = inputs.distr_params_loc_res
    u[:distr_fn_loc_work] = get(opts, inputs.distr_fn_loc_work, distr_zero)
    u[:distr_params_loc_work] = inputs.distr_params_loc_work
    u[:distr_fn_nonloc] = get(opts, inputs.distr_fn_nonloc, distr_zero)
    u[:distr_params_nonloc] = inputs.distr_params_nonloc

    println("using contact functions ", u[:distr_fn_hh], " ", u[:distr_fn_non_hh], " ", u[:distr_fn_loc_res], " ", u[:distr_fn_loc_work], " ", u[:distr_fn_nonloc])
    println("using contact params ", u[:distr_params_hh], " ", u[:distr_params_non_hh], " ", u[:distr_params_loc_res], " ", u[:distr_params_loc_work], " ", u[:distr_params_nonloc], " ", inputs.flags)

    u[:h_test] = in(:h_test,inputs.flags)
    u[:w_test] = in(:w_test,inputs.flags)

    ## queue initial events
    if haskey(inputs.init_inf, id)
        for i in rand(inputs.agents_assigned[id], inputs.init_inf[id])
            q_event!(u, infectionEvent(1,i))
        end
    end
    q_event!(u, periodicStuff(u[:periodic_stuff_period]))
    q_event!(u, reportingEvent(1))

    return nothing
end


