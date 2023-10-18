##
## note, everything here is algorithm framework, no simulation logic allowed
##

using Distributed

## the definitions in this block are sent to all processors
@everywhere begin

using DataStructures ## provides PriorityQueue

## define a simple type hierarchy for handing different kinds of events
abstract type simEvent end

abstract type abstractSimData end ## concrete type must be defined in sim logic

mutable struct SimUnit##{Ci<:Channel,Cr<:Channel}
    ## these entries are required by the worker loop
    id::Int
    inchan::RemoteChannel##{Ci}
    outchan::RemoteChannel##{Ci}
    reportchan::RemoteChannel##{Cr}
    t_inc::UInt32 #sync delay; must be specified (e.g., equal to duration of "exposed" state)
    q::PriorityQueue{simEvent,UInt32} ## event queue
    global_events::Vector{simEvent} ## global event cache
    #log::Vector{Tuple{UInt32,Any}}
    ## a struct for all other data; must be defined and initialized in sim logic
    d::abstractSimData
    ## constructor
    ##function SimUnit(id::Int,inchan::RemoteChannel{Ci},outchan::RemoteChannel{Ci},reportchan::RemoteChannel{Cr}) where {Ci<:Channel,Cr<:Channel}
    ##    return new{Ci,Cr}(id,inchan,outchan,reportchan,UInt32(5),PriorityQueue{simEvent,UInt32}(),simEvent[],Tuple{UInt32,Any}[])
    ##end
    function SimUnit(id::Int,inchan::RemoteChannel,outchan::RemoteChannel,reportchan::RemoteChannel)
        return new(id,inchan,outchan,reportchan,UInt32(5),PriorityQueue{simEvent,UInt32}(),simEvent[])
    end
end

## events should have field t for time of event, but if you write one that doesn't, just override this fn
function timeof(e::T)::UInt32 where T<:simEvent
    return e.t
end

function q_event!(u::SimUnit, e::T) where T<:simEvent
    enqueue!(u.q, e, timeof(e))
end

## special event for triggering global sync
## (this is the simplest way to avoid dequeueing events that are scheduled after events from global sync)
mutable struct syncEvent <: simEvent
    const t::UInt32
end

## special handler for sync event; performs sync and creates next sync event
## take!() will pause the worker loop until inchan returns something
function handle_event!(local_sim::SimUnit, e::syncEvent)::Vector{simEvent}
    outchan::RemoteChannel = local_sim.outchan
    inchan::RemoteChannel = local_sim.inchan
    ## note, currently put! blocks if global syncer hasn't taken the last thing we sent to outchan (shouldn't happen)
    ## sending global events only on sync, as a collection
    ##  (otherwise how will the global syncer know when all workers have processed a time period?)
    ## if outchan's memory is not on this worker, a copy is automatically made
    ## put!(outchan, local_sim.global_events)
    ## if outchan is local to this worker, need to make a copy:
    put!(outchan, deepcopy(local_sim.global_events))
    ## either way, ok to delete:
    empty!(local_sim.global_events) 
    ##await incoming global events
    events::Vector{simEvent} = take!(inchan)
    println("t = ", timeof(e), "; q len = ", length(local_sim.q), "; took ", length(events), " events from global")
    ##queue the next sync
    q_event!(local_sim, syncEvent(timeof(e) + local_sim.t_inc))
    ## distributed garbage collection is not very smart
    GC.gc()
    ##return global events
    return events
end

## special method just for events returned by sync; always place them in the local queue
function sort_events!(local_sim::SimUnit, orig::syncEvent, events::Vector{simEvent})
    for e in events
        q_event!(local_sim, e)
    end
    return nothing
end

## this is the main worker loop that runs local sim processes
## spawning this function is also how sim data (inside "unit") is sent to a processor
## (refs to inter-process comm channels are also inside "unit")
function process!(unit::SimUnit, tStop::UInt32)

    ## queue the first sync event (next will be queued by event handler)
    q_event!(unit, syncEvent(unit.t_inc))
    t = UInt32(0)

    while t < tStop
        ## get the next event in the queue, advance time
        ## note, the queue will never be empty because next sync event is always added
        e, t = dequeue_pair!(unit.q)	
        ## handle event; can change state, should generate and return future event(s)
        future_events = handle_event!(unit, e)
        ## then, either add future event(s) to local queue or save to global sync cache
        sort_events!(unit, e, future_events)
        ## note, no data reporting directly in this loop (it's an event loop, not a time loop)
    end

    println("done")
    ## close channels here so global syncer will know this worker is done
    close(unit.inchan); close(unit.outchan); close(unit.reportchan)
    return summarize(unit) ## summarize() must be defined with sim logic
end

end ## @everywhere begin


abstract type abstractGlobData end ## concrete type must be defined with sim logic

## receives global events from remote channel, returns a dict specifying which sim unit gets which events
function process_glob_channel(glob_data::U, targets::Vector{Int64}, ch::RemoteChannel{V})::Dict{Int64,Vector{simEvent}} where U<:abstractGlobData where V<:Any
    try
        global_events::Vector{simEvent} = take!(ch)
        return sort_glob_events(glob_data, targets, global_events) ## must be defined along with event handlers
    catch e
        return Dict{Int64,Vector{simEvent}}() ## empty dict if channel closed
    end
end

## create sim units, initialize with data and comm channels, and send to processors
##   tStop: remote workers need to know when to stop
##   puts handles to remote processes in "futures"
function spawn_units!(tStop::UInt32, unit_ids::Vector{Int64}, futures::Dict{Int64, Future}, in_chans::Dict{Int64, RemoteChannel{U}}, out_chans::Dict{Int64, RemoteChannel{U}}, report_chans::Dict{Int64, RemoteChannel{V}}; kwargs...) where {U,V}

    inputs = modelInputs(unit_ids; kwargs...)

    for id in unit_ids
        u = SimUnit(id, in_chans[id], out_chans[id], report_chans[id])
        init_sim_unit!(u, inputs) ## this fn adds domain-specific data and queues initial events
        ## start process; this sends data to remote worker (note, data should not remain stored in memory of main process)
        futures[id] = remotecall(process!, id, u, tStop)
    end

    return globalData(inputs) ## return data needed to sort global events among sim units
end

## "main" loop that handles global sync
## do not put any sim logic in there
function run(tStop::Int, unit_ids::Vector{Int64}; kwargs...)

    ## channels used for global sync
    ## by default, these channels' memory exists on the local process; when a remote process put!s to the
    ##  channel, the data is copied from remote memory to local memory. When the local process
    ##  put!s to a local channel, nothing is copied until the remote process take!s it. Therefore,
    ##  when sending a collection to a local channel, you must make a copy if you plan to alter the collection 
    ## currently, making all chanels exist remotely using extra argument to constructor:
    in_chans = Dict(i => RemoteChannel(()->Channel{Vector{simEvent}}(1), i) for i in unit_ids)
    out_chans = Dict(i => RemoteChannel(()->Channel{Vector{simEvent}}(1), i) for i in unit_ids)
    ## don't use report chan unless needed; more memory copying
    report_chans = Dict(i => RemoteChannel(()->Channel{Any}(100), i) for i in unit_ids)
    ## spawning a remote process returns a Future; store them here
    futures = Dict{Int64, Future}() 

    glob_data = spawn_units!(UInt32(tStop), unit_ids, futures, in_chans, out_chans, report_chans; kwargs...)
    GC.gc()
    
    ## handle global sync
    glob_dicts = Vector{Dict{Int64,Vector{simEvent}}}(undef,length(unit_ids))
    global_events = Dict{Int64,Vector{simEvent}}()
    while true
        ## sort and store global events coming from each worker
        ## asyncmap() conveniently does this as each channel reports in
        asyncmap!(i->process_glob_channel(glob_data, unit_ids, out_chans[i]), glob_dicts, unit_ids)
        ## when all workers have reported, merge results into a single dict
        mergewith!(vcat, global_events, glob_dicts...)
        ## empty dict means all channels were closed (vs no events would produce a dict of empty vectors)
        if isempty(global_events)
            break
        end

        ## send to designated workers (no need to do this async)
        for (k::Int64, v::Vector{simEvent}) in global_events
            #println("putting ", length(v), " events on inchan ", k)
            try
                ## if channel is local, probably need to make a copy
                ## or maybe not, because we're not doing anything to "v" but only to the Dict its reference lives in?
                ## put!(in_chans[k], deepcopy(v))
                ## if channel is on remote worker, data is copied automatically:
                put!(in_chans[k], v)
            catch e
                println("channel ", k, " prematurely closed")
            end
        end
        ## using mergewith!(), so empty for next loop:
        empty!(global_events)
        ## ugh distributed GC
        GC.gc()
    end
    
    ## futures contain data returned by process!()
    res = Dict{Int64,Any}()
    for (k,v) in futures
        wait(v)
        res[k] = fetch(v)
        finalize(v)
    end

    for i in unit_ids
        ## not sure it's actually necesary to "finalize" channels but doesn't seem to hurt
        remotecall_fetch(finalize, i, in_chans[i])
        remotecall_fetch(finalize, i, out_chans[i])
        remotecall_fetch(finalize, i, report_chans[i])
        #finalize(in_chans[i])
        #finalize(out_chans[i])
        #finalize(report_chans[i])
        remotecall_fetch(GC.gc, i) ## distributed GC is _really_ not very smart
    end
    GC.gc()
    return res
end

