##
## single-process version
##

include("utils.jl")
include("fileutils.jl")

using LinearAlgebra
using Random
using Distributions
using DataStructures ## provides PriorityQueue

## samples with replacement from collection a
## if a is empty, returns an empty vector instead of an error
function rsamp(a::AbstractArray{T}, n::I)::Vector{T} where {T<:Any,I<:Integer}
    return isempty(a) ? T[] : rand(a,n)
end

## struct that stores data within a SimUnit
## (this could have just been a single dict, but this way the types are identified, helps the compiler)
mutable struct SimData
    netw::Dict{Symbol,SparseMatrixCSC{Bool,UInt32}} ## together in a dict so they can be treated as a collection
    sets::Dict{Symbol,Set{UInt32}} ## which agents belong to this sim unit, along with any subsets
    pInf0::Float64 ## baseline (mean) infection probability
    pInfAmp::Float64 ## amplitude of seasonal change in infection prob
    w_pInf::Float64 ## period
    phi_pInf::Float64 ## phase
    inf_probs::Dict{Symbol,Float64} ## infection probability multipliers, stored in a dict so a different one can be associated with each network
    distr_fns::Dict{Symbol,Function} ## ditto for fns used to determine contact events
    distr_params::Dict{Symbol,Tuple} ## params for those fns
    loc_lookup::Dict{Symbol,Dict{UInt32,UInt32}} ## for indexing ephemeral contact networks by location

    n_agents::Int
    I_set::Set{UInt32} ## current infected
    R_set::Set{UInt32} ## current recovered
    t_inc::Int ## incubation time
    t_recovery::UnitRange{Int} ## range of recovery times
    pResist::Float64 ## pop-wide p prior immunity (1 - S0/N)
    report_freq::UInt32
    periodic_stuff_period::UInt32 ## mainly for updating out-of-network infections

    vaccinated::Dict{Int,Set{UInt32}} ## by _week_ : set of vaccinated individuals
    last_vacc_week::Int
    VE_inf::Float64 ## vacc effectiveness vs infection
    VE_hosp::Float64 ## vacc effectiveness vs hospitalization
    pHosp::Vector{Float64} ## p hospitalization, by group

    dummies_assigned::Vector{UInt32} ## workers without households; special logic needed; vector because randsubseq can't handle Set
    outw_assigned::Vector{UInt32} ## workers without workplaces
    mean_hh_connections::Float64 ## mainly for the dummies
    mean_wp_connections::Float64 ## mainly for out-workers

    report_groups::Set{Symbol} ## sets for which data should be reported
    cumI::Dict{Symbol,Int} ## cumulative infection counts for each set
    cumHosp::Dict{Symbol,Int} ## " hospitalizations
    report_series::Dict{Symbol,Vector{Int}} ## dict of time series for reporting
    #geo_assigned::Dict{UInt32,UInt16} ## for tracking # infections by home location
    #cumI_by_geo::Dict{Symbol,Vector{Int}} ## each location is a vector index
    #report_by_geo::Dict{Symbol,Vector{Vector{Int}}} ## time series by geo; each location is a vector index

    intervals::Dict{Symbol, UnitRange{Int}} ## time intervals for special treatment; holidays etc
    flags::Set{Symbol} ## misc boolean flags, mainly for testing
    ## constructor returns an empty struct; filled using external inputs by init_sim_unit fn below
    SimData() = new()
end


## define a simple type hierarchy for handing different kinds of events
abstract type simEvent end

## simunit consists of event queue and data
mutable struct SimUnit
    q::PriorityQueue{simEvent,UInt32} ## event queue
    d::SimData
end

## events should have field t for time of event, but if you write one that doesn't, just override this fn
function timeof(e::T)::UInt32 where T<:simEvent
    return e.t
end

function q_event!(u::SimUnit, e::T) where T<:simEvent
    enqueue!(u.q, e, timeof(e))
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
    const isVacc::Bool ## vaccinated prior to this event?
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


##
## event handlers, dispatch on event type
## worker loop expects these functions to return future events as a Vector{simEvent}
##

## most efficient way to do logging/reporting:
## add code to existing events (e.g., increment case count in E->I event)
##  and periodically queue a "data summary" event

## called when an agent enters I compartment to update cumI counts 
## if agent was also hospitalized, updates cumHosp
##   keep tracks of counts for each population subgroup in report_groups
function update_counts!(d::SimData, i::UInt32, hosp::Bool)
    ## agent's home location (= 0 for dummies, who don't have a home location)
    #g = d.geo_assigned[i]
    for k in d.report_groups
        if i in d.sets[k] ## agent id belongs to group k
            d.cumI[k] += 1 ## update pop-wide count for group k
            hosp && (d.cumHosp[k] += 1)
            #(g > 0) && (d.cumI_by_geo[k][g] += 1) ## update geo vector for k
        end
    end
end

##
## reporting event; writes data to report_series at specified time intervals
##
function handle_event!(u::SimUnit, e::reportingEvent)::Vector{simEvent}
    println("t = ",e.t," q len =",length(u.q))
    d = u.d
    ## append current count values to report series
    push!(d.report_series[:active], length(d.I_set))
    push!(d.report_series[:cumI], d.cumI[:agents_assigned])
    push!(d.report_series[:cumHosp], d.cumHosp[:agents_assigned])
    push!(d.report_series[:cumI_65o], d.cumI[:age_65o])
    push!(d.report_series[:cumHosp_65o], d.cumHosp[:age_65o])
    ## to report vectors, copy the vector to capture current values
    #for k in d.report_groups
    #    push!(d.report_by_geo[k], copy(d.cumI_by_geo[k])) 
    #end
    ## queue for next period
    q_event!(u, reportingEvent(e.t + d.report_freq))
    return simEvent[]
end

function simplecounts(v::AbstractArray{<:Integer}, n::Integer)
    r = zeros(Int,n)
    for i in v
        (0 < i <= n) && (r[i] += 1)
    end
    return r
end

##
## tells local process what to return at the end (returning the whole simunit might take too much memory)
##
function summarize(u::SimUnit)
    d = u.d
    #max_geo_index = lastindex(first(values(d.cumI_by_geo)))
    return merge(
            ## all report series
            d.report_series,
            ## geo data for each reporting group
            #Dict(Symbol("geo_"*string(k)) => d.report_by_geo[k] for k in d.report_groups),
            ## number of agents (in this unit) for each reporting group
            Dict(Symbol("n_"*string(k)) => length(d.sets[k]) for k in d.report_groups),
            ## and also by geo
            #Dict(Symbol("n_geo_"*string(k)) => simplecounts([d.geo_assigned[i] for i in d.sets[k]], max_geo_index) for k in d.report_groups),
            ## other stuff
            Dict(:endR => length(d.R_set),
                :q_len => length(u.q),
                :S0 => d.pResist,
                :VE_inf => d.VE_inf,
                :pHosp => d.pHosp
                ))
end


function infected_or_recovered(i::UInt32, e::simEvent, d::SimData)::Bool
    return in(i, d.I_set) || in(i, d.R_set) ## testing set membership is O(1)
end

## prior immunity determined randomly upon first exposure event
function prior_immunity(i::UInt32, e::simEvent, d::SimData)::Bool
    return rand() < d.pResist
end

function will_be_hosp(i::UInt32, d::SimData, isVacc::Bool)::Bool
    risk_cat = in(i, d.sets[:age_65o]) ? 2 : 1
    ### (probability of hospitalization) and not (vaccine effective)
    return (rand() < d.pHosp[risk_cat]) && !(isVacc && (rand() < d.VE_hosp))
end

## agent becomes infected (exposed)
## return an event for the target to become contagious at the correct time
## (note if an agent becomes infected twice, the earlier one will take effect & the later will do nothing)
function handle_event!(u::SimUnit, e::infectionEvent)::Vector{simEvent}
    d = u.d
    i = e.agentid
    if infected_or_recovered(i, e, d)
        return simEvent[] ## agent not susceptible
    elseif prior_immunity(i, e, d)
        push!(d.R_set, i) ## immediately add to recovered
        return simEvent[] ## no future events
    else
        ## check if patient vaccinated by time of exposure 
        ## then pass the info to becomeContagious so consequences are processed at the right time
        ##
        ## TODO: currently ignoring seroconversion time
        ##
        wk = min(1+floor(Int,e.t/7), d.last_vacc_week)
        isVacc = in(i, d.vaccinated[wk])
        return [becomeContagious(e.t + d.t_inc, i, isVacc)]
    end
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
function distr_const(neigh::AbstractVector{I}, params::Tuple{Int})::Vector{UInt32} where I<:Integer
    return rsamp(neigh, params[1])
end

## Poission(L) gives the number of contact events, then sample with replacement; L is the first/only value in params
function distr_pois(neigh::AbstractVector{I}, params::Tuple{R})::Vector{UInt32} where {I<:Integer,R<:Real}
    n = rand(Poisson(params[1]))
    return rsamp(neigh, n)
end

## Geometric (discrete analogue of Exponential) gives the number of contact events; params[1] is the mean
function safe_rgeo(u::R) where R<:Real
    return u > 0 ? rand(Geometric(1/(1+u))) : 0 ## must be > 0
end
function distr_exp(neigh::AbstractVector{I}, params::Tuple{R})::Vector{UInt32} where {I<:Integer,R<:Real}
    n = safe_rgeo(params[1])
    return rsamp(neigh, n)
end

## fn to indicate no contacts of a certain type
function distr_zero(neigh::AbstractVector{I}, params::T)::Vector{UInt32} where {I<:Integer, T<:Tuple}
    return UInt32[]
end

## look up agent's home/work location for local ephemeral contacts
function loc_res(d::SimData, i::UInt32)::UInt32
    return get(d.loc_lookup[:res], i, UInt32(0))
end
function loc_work(d::SimData, i::UInt32)::UInt32
    return get(d.loc_lookup[:work], i, UInt32(0))
end

##
## uses distr_fn to determine which network contacts become infected
##
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
## change state, and generate all future infection events that result
##
function handle_event!(u::SimUnit, e::becomeContagious)::Vector{simEvent}
    d = u.d
    i = e.agentid
    ## agent avoids illness (and does not transmit) if vaccine effective
    if (e.isVacc && (rand() < d.VE_inf))
        ## assume that "vaccine effectiveness" refers to future infection attempts also
        push!(d.R_set, i) ## by immediately marking as "recovered"
        return simEvent[] ## no future events
    else
        push!(d.I_set, i) ## append to current infected set
        hospitalized = will_be_hosp(i, d, e.isVacc) ## hospitalization currently does nothing except update the count
        update_counts!(d, i, hospitalized)
        duration = rand(d.t_recovery)
        recov_event = becomeRecovered(e.t + duration, i)
        ## collect references to networks defined in the simunit (we will broadcast the infection fn over these)
        networks = [d.netw[k] for k in [:hh, :wp, :sch, :gq, :loc_matrix_res, :loc_matrix_work]]
        ## column = agent id, except location networks where it's a location index, and nonlocal "network" that's just 1 column
        col_idxs = [i, i, i, i, loc_res(d,i), loc_work(d,i)]
        ## distribution function for # of contacts in each network; currently using non-hh for gq's
        distr_fns = [d.distr_fns[k] for k in [:hh, :non_hh, :non_hh, :non_hh, :loc_res, :loc_work]]
        distr_params = [d.distr_params[k] for k in [:hh, :non_hh, :non_hh, :non_hh, :loc_res, :loc_work]]  

        ## infection p = pInf0 + pInfAmp * cos(t*w_pInf + phi_pInf)
        pInf = d.pInf0 + d.pInfAmp * cos(e.t * d.w_pInf + d.phi_pInf)
        ## adjust by network type
        inf_ps::Vector{Float64} = pInf .* [d.inf_probs[k] for k in [:p_inf_hh, :p_inf_wp, :p_inf_sch, :p_inf_gq, :p_inf_loc, :p_inf_loc]]

        ## make modifications based on whichever special time intervals are currently active
        active_keys::Vector{Symbol} = [k for (k, i_range) in d.intervals if in(e.t, i_range)]
        for i_key in active_keys
            ## switch on i_key
            if (i_key == :sch_closed)
                ## schools are closed; disable school network for students and work networks for teachers
                col_idxs[3] = 0 ## school network only exists for k12 students
                in(i, d.sets[:k12_workers]) && (col_idxs[[2,6]] .= 0)
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

## recovery event just changes state; currently no future event
function handle_event!(u::SimUnit, e::becomeRecovered)::Vector{simEvent}
    delete!(u.d.I_set, e.agentid) ## deleting from a set is O(1)
    push!(u.d.R_set, e.agentid) ## add to recovered set
    return simEvent[]
end

## periodStuff currently just handles out-of-network infections 
##    (people who work or live outside of the synth area)
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

    ## infection p = pInf0 + pInfAmp * cos(t*w_pInf + phi_pInf)
    pInf = d.pInf0 + d.pInfAmp * cos(e.t * d.w_pInf + d.phi_pInf)
    ## adjust by netw type
    p_not_w = 1.0 - pInf * d.inf_probs[:p_inf_wp]
    p_not_h = 1.0 - pInf * d.inf_probs[:p_inf_hh]
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

## no sorting in this version, just add to queue
function sort_events!(local_sim::SimUnit, orig::simEvent, events::Vector{simEvent})
    for e in events
        q_event!(local_sim, e)
    end
    return nothing
end


## mean household connections, excluding household-less dummies
calc_mean_hh(netw_mat,exclude) = mean(sum(netw_mat[:,setdiff(axes(netw_mat,2), exclude)], dims=1))
## mean workplace connections = mean non-hh cnxs for people with 1+ non-hh cnxs, excluding workplace-less out-workers
calc_mean_wp(netw_mat,exclude) = mean(filter(x->x>0, sum(netw_mat[:,setdiff(axes(netw_mat,2), exclude)], dims=1)))


function modelInputs(;kwargs...)
    ## note, calling kwargs like a fn is enabled by a definition in utils.jl
    netw_hh = kwargs(:netw_hh, ()->SparseMatrixCSC{Bool, UInt32}(Symmetric(dser_path("jlse/adj_mat_hh.jlse"))))
    netw_wp = kwargs(:netw_wp, ()->SparseMatrixCSC{Bool, UInt32}(Symmetric(dser_path("jlse/adj_mat_wp.jlse"))))
    ## dummies appear in work networks, but live outside the synth pop (no household or demographic info generated)
    ## local sim unit is responsible for determining if they got infected at home
    ## collect() because randsubseq can't handle sets
    dummies = kwargs(:dummies, ()->collect(UInt32, keys(dser_path("jlse/adj_dummy_keys.jlse"))))
    ## outside workers have households, but no workplace network
    ## local sim unit is responsible for determining if they got infected at work
    outw = kwargs(:out_workers, ()->collect(UInt32, keys(dser_path("jlse/adj_out_workers.jlse"))))
    
    println("returning model inputs")
    return Dict(
        :netw_hh => netw_hh,
        ## must be something wih O(1) lookup:
        :agents_assigned => kwargs(:agents_assigned ,()->Set(UInt32.(collect(1:size(netw_hh,2))))),
        :netw_wp => netw_wp,
        :netw_sch => kwargs(:netw_sch, ()->SparseMatrixCSC{Bool, UInt32}(Symmetric(dser_path("jlse/adj_mat_sch.jlse")))),
        :netw_gq => kwargs(:netw_gq, ()->SparseMatrixCSC{Bool, UInt32}(Symmetric(dser_path("jlse/adj_mat_gq.jlse")))),
        :loc_matrix_res => kwargs(:loc_matrix_res, ()->SparseMatrixCSC{Bool,UInt32}(dser_path("jlse/res_loc_contact_mat.jlse"))),
        :loc_matrix_work => kwargs(:loc_matrix_work, ()->SparseMatrixCSC{Bool,UInt32}(dser_path("jlse/work_loc_contact_mat.jlse"))),
        :loc_lookup_res => kwargs(:loc_lookup_res, ()->Dict{UInt32,UInt32}(dser_path("jlse/res_loc_lookup.jlse"))),
        :loc_lookup_work => kwargs(:loc_lookup_work, ()->Dict{UInt32,UInt32}(dser_path("jlse/work_loc_lookup.jlse"))),
        :geo_lookup => kwargs(:geo_lookup, ()->Dict{UInt32,UInt16}(dser_path("precalc/cbg_idx_lookup.jlse"))),
        :t_inc => UInt32(get(kwargs, :t_inc, 2)),
        :t_recovery => get(kwargs, :t_recovery, 4:8),
        :pResist => get(kwargs, :pResist, 0.333),
        :init_inf => get(kwargs, :init_inf, 100),
        :dummies => dummies,
        :out_workers => outw,
        :k12_workers => kwargs(:k12_workers, ()->Set{UInt32}(dser_path("precalc/k12_workers.jlse"))),
        :age_65o => kwargs(:age_65o, ()->Set{UInt32}(dser_path("precalc/age_65o.jlse"))),
        :vaccinated => kwargs(:vaccinated, ()->dser_path("precalc/vaccinated.jlse")),
        :VE_inf => get(kwargs, :VE_inf, 0.25),
        :pHosp => get(kwargs, :pHosp, [0.01, 0.1]),
        :mean_hh_connections => kwargs(:mean_hh_connections, ()->calc_mean_hh(netw_hh,dummies)),
        :mean_wp_connections => kwargs(:mean_wp_connections, ()->calc_mean_wp(netw_wp,outw)),
        :report_freq => UInt32(get(kwargs, :report_freq, 7)),
        ## seasonal infection probability; assuming annual period
        :pInf0 => get(kwargs, :pInf0, 0.1), ## baseline(mean)
        :pInfAmp => get(kwargs, :pInfAmp, 0.05), ## amplitude
        :tPeak => get(kwargs, :tPeak, 120), ## days from sim start to peak inf prob
        ## p_inf's below are multipliers
        ## defaults for within-household infection
        :p_inf_hh => get(kwargs, :p_inf_hh, 1.0), 
        :distr_fn_hh => get(kwargs, :distr_fn_hh, :all),
        :distr_params_hh => get(kwargs, :distr_params_hh, ()),
        ## defaults for work,school,GQ infection
        :p_inf_wp => get(kwargs, :p_inf_wp, 1.0),
        :p_inf_sch => get(kwargs, :p_inf_sch, 1.0),
        :p_inf_gq => get(kwargs, :p_inf_gq, 1.0),
        :distr_fn_non_hh => get(kwargs, :distr_fn_non_hh, :all),
        :distr_params_non_hh => get(kwargs, :distr_params_non_hh, ()),
        ## defaults for ephemeral/location-based infection
        :p_inf_loc => get(kwargs, :p_inf_loc, 1.0),
        :distr_fn_loc_res => get(kwargs, :distr_fn_loc_res, :zero),
        :distr_params_loc_res => get(kwargs, :distr_params_loc_res, ()),
        :distr_fn_loc_work => get(kwargs, :distr_fn_loc_work, :zero),
        :distr_params_loc_work => get(kwargs, :distr_params_loc_work, ()),
        :distr_fn_nonloc => get(kwargs, :distr_fn_nonloc, :zero),
        :distr_params_nonloc => get(kwargs, :distr_params_nonloc, ()),
        ## special time intervals; missing by default
        :sch_closed => get(kwargs, :sch_closed, missing),
        ## misc flags
        :flags => kwargs(:flags, ()->Set(Symbol[]))
        )
end

## queue initial infections when init_inf is # infections
function q_init_inf!(u::SimUnit, init_inf::I) where I<:Integer
    for i in rand(collect(u.d.sets[:agents_assigned]), init_inf)
        #q_event!(u, infectionEvent(1,i))
        q_event!(u, becomeContagious(0,i,false))
    end
end

## queue initial infections when init_inf is a list of agent ids
function q_init_inf!(u::SimUnit, init_inf::Vector{I}) where I<:Integer
    for i in init_inf
        if i in u.d.sets[:agents_assigned]
            #q_event!(u, infectionEvent(1,i))
            q_event!(u, becomeContagious(0,i,false))
        end
    end
end


## simunit constructor, from modelinputs dict
## adds initial event(s) to queue (at least one unit must have an initial event)
function SimUnit(inputs::Dict{Symbol,Any})

    q = PriorityQueue{simEvent,UInt32}()
    d = SimData()

    d.n_agents = length(inputs[:agents_assigned])
    d.I_set = Set{UInt32}()
    d.R_set = Set{UInt32}()
    d.t_inc = inputs[:t_inc] 
    d.t_recovery = inputs[:t_recovery]
    ## infection probability is defined in terms of infectiousness duration, so 
    ##  update out-of-network infections on the same timescale so probabilities are correct
    d.periodic_stuff_period = round(UInt32, mean(inputs[:t_recovery]))
    d.report_freq = inputs[:report_freq]

    ## prior immunity is decided when needed
    d.pResist = inputs[:pResist]

    ## seasonal infection probability; assuming annual period
    ## p = pInf0 + pInfAmp * cos(t*w_pInf + phi_pInf)
    d.pInf0 = inputs[:pInf0] ## baseline(mean)
    d.pInfAmp = inputs[:pInfAmp] ## amplitude
    tPeak = inputs[:tPeak] ## days from sim start to peak inf prob
    d.w_pInf = 2*pi/365.24
    d.phi_pInf = -2*pi*tPeak/365.24

    d.vaccinated = inputs[:vaccinated] ## by _week_ : set of vaccinated individuals
    d.last_vacc_week = maximum(keys(d.vaccinated))
    d.VE_inf = inputs[:VE_inf] ## vacc effectiveness vs infection
    d.VE_hosp = 0.5 ## vacc effectiveness vs hospitalization; currently fixed at 50%
    d.pHosp = inputs[:pHosp] ## p hospitalization, by group; currently ages 64- and 65+

    ## infection probability multipliers for each type of network
    d.inf_probs = Dict(
        :p_inf_hh => inputs[:p_inf_hh],
        :p_inf_wp => inputs[:p_inf_wp],
        :p_inf_sch => inputs[:p_inf_sch],
        :p_inf_gq => inputs[:p_inf_gq],
        :p_inf_loc => inputs[:p_inf_loc]
    )

    ## having many network matrices wastes memory because indices are repeated; could combine into an 8-bit binary
    d.netw = Dict(
        :hh => inputs[:netw_hh],
        :wp => inputs[:netw_wp],
        :sch => inputs[:netw_sch],
        :gq => inputs[:netw_gq],
        :loc_matrix_res => inputs[:loc_matrix_res], ## local ephemeral contacts by location; every sim unit needs the full matrix
        :loc_matrix_work => inputs[:loc_matrix_work],
        ## non-local contact network is everyone in the population, so just a single column of "trues" 
        ## TODO: this wastes memory; but we're not currently using it
        ## :nonlocal => SparseMatrixCSC{Bool, UInt32}(trues(size(inputs[:netw_hh], 1), 1))
    )

    ## agent assignment and other membership groups
    ##   if there are many of these, could change representation to bitmatrix or something 
    d.sets = merge(
        Dict(:agents_assigned => inputs[:agents_assigned]), 
        Dict(k => inputs[k] for k in 
            [:k12_workers,:age_65o]))

    ## keys from d.sets for which data should be collected (:agents_assigned is everyone in this sim unit)
    d.report_groups = Set([:agents_assigned,
                            :age_65o])
    ## keep track of cumulative infections for those groups
    d.cumI = Dict(d.report_groups .=> 0)
    d.cumHosp = Dict(d.report_groups .=> 0)
    ## initialize time series to be updated in reporting event handler
    d.report_series = Dict(k => Int[] for k in [:active, :cumI, 
                                                :cumHosp, :cumI_65o, :cumHosp_65o])

    ## for tracking # infections by home location
    #d.geo_assigned = inputs[:geo_lookup]
    #max_geo_index = maximum(values(inputs[:geo_lookup]))
    ## do this separately for each reporting group
    #d.cumI_by_geo = Dict(k=>zeros(Int, max_geo_index) for k in d.report_groups)
    #d.report_by_geo = Dict(k=>Vector{Int}[] for k in d.report_groups) ## a vector of vectors: [timepoint][geo]

    ## note, loc_lookup_res is missing inst gq residents (and dummies); loc_lookup_work has only commuters
    d.loc_lookup = Dict(
        :res => inputs[:loc_lookup_res],
        :work => inputs[:loc_lookup_work]
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

    d.dummies_assigned = inputs[:dummies]
    d.outw_assigned = inputs[:out_workers]
    ## note, these are currently global means (not per sim unit)
    d.mean_hh_connections = inputs[:mean_hh_connections]
    d.mean_wp_connections = inputs[:mean_wp_connections]

    ## include keys only for specified intervals
    d.intervals = Dict{Symbol, UnitRange{Int}}(k=>inputs[k] for k in 
        [:sch_closed] if !ismissing(inputs[k]))
    
    d.flags = inputs[:flags]

    u = SimUnit(q,d)
    ## queue initial events
    q_init_inf!(u, inputs[:init_inf])
    q_event!(u, periodicStuff(d.periodic_stuff_period))
    q_event!(u, reportingEvent(d.report_freq))

    println("using contact functions ", d.distr_fns[:hh], " ", d.distr_fns[:non_hh], " ", d.distr_fns[:loc_res], " ", d.distr_fns[:loc_work])
    println("using contact params ", d.distr_params[:hh],d.inf_probs[:p_inf_hh], " ", d.distr_params[:non_hh],d.inf_probs[:p_inf_wp],d.inf_probs[:p_inf_sch],d.inf_probs[:p_inf_gq], " ", d.distr_params[:loc_res],d.inf_probs[:p_inf_loc], " ", d.distr_params[:loc_work],d.inf_probs[:p_inf_loc], " ", d.flags)
    return u
end


## event loop
function run(unit::SimUnit, tStop::I) where I<:Integer

    t = UInt32(0)
    while t < tStop
        ## get the next event in the queue, advance time
        ## note, the queue should never be empty
        e, t = dequeue_pair!(unit.q)
        ## handle event; can change state, should generate and return future event(s)
        future_events = handle_event!(unit, e)
        sort_events!(unit, e, future_events)
        ## note, no data reporting directly in this loop (it's an event loop, not a time loop)
    end

    return summarize(unit) ## summarize() must be defined with sim logic
end



##
## TODO: don't allow decrease (ipf doesn't know this rule)
##
##
## TODO: can make this a lot more memory-efficient for a small runtime cost 
##
function assign_vac(groups, vac_wk_loc_age, assign_outside, vac_outside_wk)
    vacSet_by_wk = Dict{Int,Set{UInt32}}()
    for wk in eachindex(vac_wk_loc_age)
        v_loc_age = vac_wk_loc_age[wk]
        v_out = vac_outside_wk[wk]
        vacSet_by_wk[wk] = Set([reduce(vcat, [
                reduce(vcat, [groups[i][loc][1:v_loc_age[loc,i]] for i in eachindex(groups)])
            for loc in axes(v_loc_age,1)]); assign_outside[1:v_out]])
    end
    return vacSet_by_wk
end


## pre-calculate stuff used in intialization


## can store people in specified industries, schools, etc., as a set of network indices
##   if using several of these sets, could save memory (with only a small performace cost) by combining
##   them with agents_assigned into a single dict
function precalc_sets()

    p_idxs = let k = dser_path("jlse/adj_mat_keys.jlse"); Dict(k .=> UInt32.(eachindex(k))); end

    mkpath("precalc")
    
    ## for looking up a person's location index
    ser_path("precalc/cbg_idx_lookup.jlse", Dict(i=>k[3] for (k,i) in p_idxs))

    ## list of person indices by residence location (keyed by geocode)
    cbg_k = dser_path("jlse/cbg_idxs.jlse") ## mapping of location index to location geocode
    cbg_k[0] = "outside"
    gdf = groupby(DataFrame([(p=k[1],hh=k[2],cbg=k[3],idx=v) for (k,v) in p_idxs]), "cbg")
    p_idxs_all_by_h_cbg = Dict(cbg_k[gk["cbg"]] => gdf[gk][:,:idx] for gk in keys(gdf))
    ser_path("precalc/p_idxs_all_by_h_cbg.jlse", p_idxs_all_by_h_cbg)

    work_in_school = let sch_workers = dser_path("jlse/sch_workers.jlse"); Set(p_idxs[k[1:3]] for k in reduce(vcat, collect(values(sch_workers)))); end

    people = dser_path("jlse/people.jlse")
    p_in_school = Set(p_idxs[k[1:3]] for k in keys(filterv(p->(!ismissing(p.sch_grade) && !in(p.sch_grade, ["c","g"])), people)))
    school_student_or_worker = union(p_in_school,work_in_school)
    age65o = Set(p_idxs[k[1:3]] for k in keys(filterv(x->x.age>=65, people)))
    age17u = Set(p_idxs[k[1:3]] for k in keys(filterv(x->x.age<18, people)))
    age_18_49 = Set(p_idxs[k[1:3]] for k in keys(filterv(x->(17<x.age<50), people)))
    age_50_64 = Set(p_idxs[k[1:3]] for k in keys(filterv(x->(49<x.age<65), people)))

    ser_path("precalc/age_17u.jlse", age17u)
    ser_path("precalc/age_18_49.jlse", age_18_49)
    ser_path("precalc/age_50_64.jlse", age_50_64)
    ser_path("precalc/age_65o.jlse", age65o)    

    ser_path("precalc/k12_students.jlse", p_in_school)
    ser_path("precalc/k12_workers.jlse", work_in_school)
    ser_path("precalc/k12_student_or_worker.jlse", school_student_or_worker)
    return nothing
end



## precalculate stuff for vaccination scenarios

using ProportionalFitting
using Dates
using Logging
using Distributions

function precalc_vacc()

    Logging.disable_logging(Logging.Info)

    p_by_cbg = dser_path("precalc/p_idxs_all_by_h_cbg.jlse")
    age_17u = dser_path("precalc/age_17u.jlse")
    age_18_49 = dser_path("precalc/age_18_49.jlse")
    age_50_64 = dser_path("precalc/age_50_64.jlse")
    age_65o = dser_path("precalc/age_65o.jlse")

    N_by_age = [length(age_17u) , length(age_18_49) , length(age_50_64) , length(age_65o)]

    d = Dict{String,Vector{UInt32}}()
    for (k,v) in p_by_cbg
        c = k[1:5]
        ppl = get(d, c, UInt32[])
        d[c] = [ppl; v]
    end
    p_by_county = Dict(k => Set(v) for (k,v) in d)
    counties = sort(collect(keys(p_by_county)))[1:end-1] ## drop "outsi"
    pop_loc = [length(p_by_county[c]) for c in counties]

    df = CSV.read("precalc/covid_vacc_by_county.csv", DataFrame)
    df.Date = Date.(df.Date,"mm/dd/yyyy");
    df = subset(df, :Date => x -> x .== Date("2023-05-10"), :FIPS => x -> x .!= "UNK", :Series_Complete_Pop_Pct => x -> .!ismissing.(x))
    df[!,:pct] = round.(df[!,:Series_Complete_Pop_Pct] / 100; digits=3)
    covid_vacc_by_fips = Dict(df.FIPS .=> df.pct)
    ser_path("precalc/covid_vacc_by_fips.jlse", covid_vacc_by_fips)
    #covid_vacc_by_fips = dser_path("precalc/covid_vacc_by_fips.jlse")
    cvac_loc = pop_loc .* [covid_vacc_by_fips[c] for c in counties]
    cvac_total = sum(cvac_loc)

    age_17u_county = [collect(intersect(p_by_county[c], age_17u)) for c in counties]
    age_18_49_county = [collect(intersect(p_by_county[c], age_18_49)) for c in counties]
    age_50_64_county = [collect(intersect(p_by_county[c], age_50_64)) for c in counties]
    age_65o_county = [collect(intersect(p_by_county[c], age_65o)) for c in counties]

    ser_path("precalc/age_17u_county.jlse",age_17u_county)
    ser_path("precalc/age_18_49_county.jlse",age_18_49_county)
    ser_path("precalc/age_50_64_county.jlse",age_50_64_county)
    ser_path("precalc/age_65o_county.jlse",age_65o_county)

    N_loc_age = [length.(age_17u_county) length.(age_18_49_county) length.(age_50_64_county) length.(age_65o_county)]

    codes = CSV.read("precalc/series_codes.csv", DataFrame, types=Dict("fips"=>String7))
    name_to_fips = Dict(codes.name .=> codes.fips)

    df = CSV.read("precalc/flu_vacc_state_age.csv", DataFrame)
    df[!,:fips] = [get(name_to_fips,k,missing) for k in df.Geography]
    df[!,:p_vacc] = df[!,"flu.coverage.rd2526.sc_A"] / 100
    df[!,:p_reduced] = df[!,"flu.coverage.rd2526.sc_B"] / 100

    vac_dates = sort(unique(df[!,:Week_Ending_Sat]))
    age_groups = ["6 Months - 17 Years","18-49 Years","50-64 Years","65+ Years"]

    df_state = subset(df, :fips=>x->mtrue.(x.=="24"))

    df_p_vacc = unstack(df_state,:Week_Ending_Sat,:Age,:p_vacc)
    ## make sure it's sorted
    sort!(df_p_vacc,[:Week_Ending_Sat])
    p_vacc_date_age = Matrix(df_p_vacc[:,age_groups])
    N_vacc_date_age = N_by_age' .* p_vacc_date_age

    df_p_reduced = unstack(df_state,:Week_Ending_Sat,:Age,:p_reduced)
    ## make sure it's sorted
    sort!(df_p_reduced,[:Week_Ending_Sat])
    p_reduced_date_age = Matrix(df_p_reduced[:,age_groups])
    N_reduced_date_age = N_by_age' .* p_reduced_date_age

    function est_vac(wk, N_vacc_date_age)
        vac_age = N_vacc_date_age[wk,:]
        fluvac_total = sum(vac_age)
        vac_loc = cvac_loc * fluvac_total / cvac_total

        M_init = vac_loc .* N_loc_age ./ pop_loc
        fac = ipf(M_init, [vac_loc, vac_age])
        M = Array(fac) .* M_init

        ## this shouldn't happen too often, but don't allow >99% vacc
        problems = findall(M .> 0.99*N_loc_age)
        (length(problems) > 0) && println("warning: wk ", wk, " vac limited by pop at ", problems)
        M = min.(M, 0.99*N_loc_age)

        return round.(Int,M)
    end

    println("estimating vac by loc and age")
    vac_wk_loc_age = map(i->est_vac(i,N_vacc_date_age), eachindex(vac_dates))
    println("estimating reduced vac by loc and age")
    vreduced_wk_loc_age = map(i->est_vac(i,N_reduced_date_age), eachindex(vac_dates))

    ser_path("precalc/p_vacc_date_age.jlse",p_vacc_date_age)
    ser_path("precalc/p_reduced_date_age.jlse",p_reduced_date_age)
    ser_path("precalc/vac_wk_loc_age.jlse",vac_wk_loc_age)
    ser_path("precalc/vreduced_wk_loc_age.jlse",vreduced_wk_loc_age)
    ser_path("precalc/p_idxs_outside.jlse", p_by_cbg["outside"])

    return nothing
end


# generate and save parameter sets
function precalc_params()
    dist_S = 0.23 + 0.15*Beta(6,3)
    dist_I = TriangularDist(500, 1500, 1000)
    dist_VE = TriangularDist(0.15, 0.35, 0.25)
    dist_H64u = TriangularDist(0.005, 0.015, 0.010)
    dist_H65o = TriangularDist(0.1, 0.3, 0.2)

    param_sets = [(S0 = rand(dist_S), I0 = round(Int, rand(dist_I)), VE = rand(dist_VE), pH64u = rand(dist_H64u), pH65o = rand(dist_H65o) )
                    for i in 1:10000];
    ser_path("param_sets.jlse", param_sets)
end



