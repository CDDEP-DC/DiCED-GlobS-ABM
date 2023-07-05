using LinearAlgebra
include("utils.jl")
include("fileutils.jl")

## remote workers need these definitions
@everywhere begin

using Random
using SparseArrays

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

function handle_event!(u::SimUnit, e::becomeContagious)::Vector{simEvent}
    ## if agent is susceptible, change state, and generate all future infection events that result
    i = e.agentid
    if in(i, u[:I_set]) || in(i, u[:R_set]) ## testing set membership is O(1)
        return simEvent[] ## agent not susceptible
    else
        u[:cum_I] += 1 ## update cumulative infection count
        push!(u[:I_set], i) ## append to current infected set
        duration = rand(u[:t_recovery])
        recov_event = becomeRecovered(e.t + duration, i)
        neigh = findall(u[:netw][:,i]) ## note, this is not really a "find", it's just a lookup in a sparse array
        if isempty(neigh)
            return [recov_event] ## no neighbors
        else
            infected = randsubseq(neigh, u[:p_inf]) ## fn provided by Random, promises efficient Bernouilli sampling
            infection_events = [infectionEvent(e.t + rand(0:duration), targ) for targ in infected]
            ## return infection events and recovery event
            return [infection_events; recov_event]
        end
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
    P_infected = length(u[:I_set]) / u[:n_agents]

    ## mean number of infected home / work connections
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

    ## append proportion infected to log
    P_infected = length(u[:I_set]) / u[:n_agents]
    push!(u[:log], (e.t,P_infected))

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
##  parametric types T U V etc can be anything the functions below are written to handle
##  (this didn't really need to be a struct, but it makes the syntax below cleaner)
@kwdef struct modelInputs{T,U,V,W} <: abstractModelInputs
    netw::T
    agents_assigned::U
    t_inc::Int64
    p_inf::Float64
    t_recovery::V
    init_inf::Dict{Int64,Int64}
    dummies_assigned::W
    outw_assigned::W
    mean_hh_connections::Float64
    mean_wp_connections::Float64
    report_freq::Int64
end

## constructor, called in gabm.jl spawn_units!()
function modelInputs(unit_ids::Vector{Int64}; kwargs...)

    ## these funcs only evaluated when needed
    read_dummies() = Set(keys(dser_path("jlse/adj_dummy_keys.jlse")))
    read_outw() = Set(keys(dser_path("jlse/adj_out_workers.jlse")))
    read_hh_net() = sparse(Symmetric(dser_path("jlse/hh_adj_mat.jlse")))
    read_non_hh() = sparse(Symmetric(dser_path("jlse/adj_mat.jlse")))
    ## mean household connections, excluding household-less dummies
    calc_mean_hh(netw_hh,dummies) = (println("calc mean hh"); mean(sum(netw_hh[:,setdiff(axes(netw_hh,2), dummies)], dims=1)))
    ## mean workplace connections = mean non-hh cnxs for people with 1+ non-hh cnxs, excluding workplace-less out-workers
    calc_mean_wp(netw_non_hh,outw) = (println("calc mean wp"); mean(filter(x->x>0, sum(netw_non_hh[:,setdiff(axes(netw_non_hh,2), outw)], dims=1))))
    ## for now, just join household and non-hh networks; keyword arg takes precedence if present
    calc_full_net(netw_hh,netw_non_hh)::SparseMatrixCSC{Bool, Int64} = (println("calc full net"); netw_hh .| netw_non_hh)

    ## awkward syntax but avoids costly evaluation when kwargs has the key
    netw_hh::SparseMatrixCSC{Bool, Int64} = get(read_hh_net, kwargs, :netw_hh)
    netw_non_hh::SparseMatrixCSC{Bool, Int64} = get(read_non_hh, kwargs, :netw_non_hh)
    ## dummies appear in work networks, but live outside the synth pop (no household or demographic info generated)
    ## local sim unit is responsible for determining if they got infected at home
    dummies::Set{Int64} = get(read_dummies, kwargs, :dummies)
    ## outside workers have households, but no workplace network
    ## local sim unit is responsible for determining if they got infected at work
    outw::Set{Int64} = get(read_outw, kwargs, :out_workers)

    mu_hh_cnx::Float64 = get(()->calc_mean_hh(netw_hh,dummies), kwargs, :mean_hh_connections)
    mu_work_cnx::Float64 = get(()->calc_mean_wp(netw_non_hh,outw), kwargs, :mean_wp_connections)

    ## assigning agents to sim units
    ## TODO: optimize
    n = length(unit_ids)
    assign_idxs = Dict(unit_ids .=> ranges(lrRound(fill(size(netw_hh,2)/n, n))))
    assign_dummies = Dict(i => collect(intersect(dummies, assign_idxs[i])) for i in unit_ids)
    assign_outw = Dict(i => collect(intersect(outw, assign_idxs[i])) for i in unit_ids)

    println("returning model inputs")
    ## this constructor syntax is enabled by @kwdef macro above
    return modelInputs(
        t_inc = get(kwargs, :t_inc, 5),
        p_inf = get(kwargs, :p_inf, 0.2),
        t_recovery = get(kwargs, :t_recovery, 8:12),
        init_inf = Dict(get(kwargs, :init_inf, [first(unit_ids) => 10])),
        netw = get(()->calc_full_net(netw_hh,netw_non_hh), kwargs, :full_net),
        agents_assigned = assign_idxs,
        dummies_assigned = assign_dummies,
        outw_assigned = assign_outw,
        mean_hh_connections = mu_hh_cnx,
        mean_wp_connections = mu_work_cnx,
        report_freq = get(kwargs, :report_freq, 5)
        )
end

## globalData constructor, called in gabm.jl spawn_units!()
## should return whatever sort_glob_events() needs
function globalData(inputs::modelInputs)
    g = zeros(UInt8,size(inputs.netw,2));
    for (k,r) in inputs.agents_assigned
        g[r] .= k
    end
    return globalData(g)
end

## this fn is called from a fn in gabm.jl
## adds domain-specific data to a sim unit
## adds initial event(s) to queue (at least one unit must have an initial event)
## (note, SimUnit type is currently just a Dict with symbol keys, can add any key needed by fns above)
function init_sim_unit!(u::SimUnit, id::Int64, inputs::modelInputs)

    println("init sim unit ", id)

    ## add state variables and parameters
    u[:netw] = spzeros(Bool,size(inputs.netw)) ## zeros use no memory in a sparse array
    u[:netw][:,inputs.agents_assigned[id]] = copy(inputs.netw[:,inputs.agents_assigned[id]]) ## copy just the columns assigned to this unit
    u[:agents_assigned] = inputs.agents_assigned[id]
    u[:n_agents] = length(u[:agents_assigned])
    u[:dummies_assigned] = inputs.dummies_assigned[id]
    u[:outw_assigned] = inputs.outw_assigned[id]
    u[:cum_I] = 0
    u[:I_set] = Set{UInt32}()
    u[:R_set] = Set{UInt32}()
    u[:t_inc] = inputs.t_inc
    u[:p_inf] = inputs.p_inf
    u[:t_recovery] = inputs.t_recovery
    ## note, these are currently global means (not per sim unit)
    u[:mean_hh_connections] = inputs.mean_hh_connections
    u[:mean_wp_connections] = inputs.mean_wp_connections
    ## infection probability is defined in terms of infectiousness duration, so 
    ##  update out-of-network infections on the same timescale so probabilities are correct
    u[:periodic_stuff_period] = round(Int, mean(inputs.t_recovery))
    u[:report_freq] = inputs.report_freq

    ## queue initial events
    if haskey(inputs.init_inf, id)
        for i in rand(inputs.agents_assigned[id], inputs.init_inf[id])
            q_event!(u, infectionEvent(1,i))
        end
    end
    q_event!(u, periodicStuff(u[:periodic_stuff_period]))
    q_event!(u, reportingEvent(u[:report_freq]))

    return nothing
end


