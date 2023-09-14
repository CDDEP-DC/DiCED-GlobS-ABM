include("gabm.jl")
include("sim_logic.jl")

function test(t; kwargs...)
    x = run(t, workers(); kwargs...)
    return Dict(k => (v[:cum_I],v[:curr_I],v[:curr_R],v[:q_len],v[:glob_len]) for (k,v) in x),
            Dict(k => v[:n_agents] for (k,v) in x),
            Dict(k => v[:log] for (k,v) in x) ## proportion infected is in here
end

#adj_mat = sparse(Symmetric(dser_path("jlse/adj_mat.jlse")));
#hh_adj_mat = sparse(Symmetric(dser_path("jlse/hh_adj_mat.jlse")));
#dum_idxs = Set(keys(dser_path("jlse/adj_dummy_keys.jlse")));
#outw_idxs = Set(keys(dser_path("jlse/adj_out_workers.jlse")));
#mu_hh_cnx = mean(sum(hh_adj_mat[:,setdiff(axes(hh_adj_mat,2), dum_idxs)], dims=1))
#mu_work_cnx = mean(filter(x->x>0, sum(adj_mat[:,setdiff(axes(adj_mat,2), outw_idxs)], dims=1)))
#netw = adj_mat .| hh_adj_mat;

t = 602
I0 = [2=>100,3=>100,4=>100]
pI = 0.12

## read network from jlse dir, run with above params, save output
for i in 2:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_constN"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
    distr_params_hh=(8,)
    ))
    GC.gc()
end

for i in 2:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2_workloc"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        distr_params_hh=(8,), distr_params_non_hh=(6,), distr_fn_loc_work=:const, distr_params_loc_work=(2,)
        ))
    GC.gc()
end

for i in 2:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2_resloc"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        distr_params_hh=(6,), distr_fn_loc_res=:const, distr_params_loc_res=(2,)
        ))
    GC.gc()
end

for i in 2:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2_res_to_glob"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        distr_params_hh=(6,), distr_fn_nonloc=:const, distr_params_nonloc=(2,)
        ))
    GC.gc()
end

for i in 2:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2_work_to_glob"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        distr_params_hh=(8,), distr_params_non_hh=(6,), distr_fn_nonloc=:const, distr_params_nonloc=(2,)
        ))
    GC.gc()
end


## generate random one-time contacts, preserving # contacts per person
#for i in 1:2
#    randonet = colshuffle(sparse(Symmetric(dser_path("jlse/adj_mat.jlse"))) .| sparse(Symmetric(dser_path("jlse/hh_adj_mat.jlse"))));
#    GC.gc()
#    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_rnet"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI, full_net=randonet))
#    GC.gc()
#end

#netw = sparse(Symmetric(dser_path("jlse/adj_mat.jlse"))) .| sparse(Symmetric(dser_path("jlse/hh_adj_mat.jlse")));
#dum_idxs = Set(keys(dser_path("jlse/adj_dummy_keys.jlse")));
#outw_idxs = Set(keys(dser_path("jlse/adj_out_workers.jlse")));
#mask = setdiff(axes(netw,2), dum_idxs, outw_idxs)
#mu = mean(sum(netw[:,mask], dims=1))
#nn = size(netw,1)
#netw, dum_idxs, outw_idxs, mask = nothing, nothing, nothing, nothing
#GC.gc()
mu = 7.793946656329233 ## MD
nn = 9440777 ## MD

## generate random one-time contacts, preserving only the mean # contacts per person
#for i in 1:2
#    really_rand = randnet(nn, mu, 1.0)
#    GC.gc()
#    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_rrand"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI, full_net=really_rand, dummies = Set(Int64[]), out_workers = Set(Int64[])))
#    GC.gc()
#end

using Graphs

## watts_strogatz "small world" network -- doesn't form hubs
for i in 1:2
    sw_net = Bool.(Graphs.LinAlg.adjacency_matrix( watts_strogatz(nn,8,0.25) ));
    GC.gc()
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_sw"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI, full_net=sw_net, dummies = Set(Int64[]), out_workers = Set(Int64[])))
    GC.gc()
end

## barabasi_albert "scale free" network -- lots of hubs
for i in 1:2
    ba_net = Bool.(Graphs.LinAlg.adjacency_matrix( barabasi_albert(nn, 4) ));
    GC.gc()
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_ba"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI, full_net=ba_net, dummies = Set(Int64[]), out_workers = Set(Int64[])))
    GC.gc()
end

## erdos_renyi "random" network
for i in 1:2
    er_net = Bool.(Graphs.LinAlg.adjacency_matrix( erdos_renyi(nn,round(Int,0.5*mu*nn)) ));
    GC.gc()
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_er"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI, full_net=er_net, dummies = Set(Int64[]), out_workers = Set(Int64[])))
    GC.gc()
end


## swap out files...
## read network from jlse dir, run with above params, save output
for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_NO_INC"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI))
    GC.gc()
end


## lower init
t = 602
I0 = [2=>10,3=>0,4=>0]
pI = 0.15
## swap out files...
## read network from jlse dir, run with above params, save output
for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_I10"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI))
    GC.gc()
end
## swap out files...
## read network from jlse dir, run with above params, save output
for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_I10_no_inc"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI))
    GC.gc()
end








using Plots

## read logs; arrange time series in matrix columns for Plots.plot()
read_log(f::String,rep::Int) = deserialize(f*string(rep)*".jlse")[3]
read_logs(f,xlim,ids,reps) = hcat([[t[2] for t in read_log(f,rep)[i][1:xlim]] for i in ids for rep in reps]...)
function read_sums(f,xlim,ids,reps)
    ## stack reps
    return hcat([ 
    sum(reduce(hcat, [[t[2] for t in read_log(f,r)[i][1:xlim]] for i in ids]),dims=2) ## sum workers
    for r in reps]...)
end

r_inc = 5
ids = 2:4

xlim = floor(Int, 600/r_inc)
xs = collect(1:xlim) .* r_inc;

ys_cN = read_sums("sim_MD_p012_constN",xlim,ids,1:3)
ys_work2 = read_sums("sim_MD_p012_2_workloc",xlim,ids,1:3)
ys_res2 = read_sums("sim_MD_p012_2_resloc",xlim,ids,1:3)
ys_w2glob = read_sums("sim_MD_p012_2_work_to_glob",xlim,ids,1:3)
ys_r2glob = read_sums("sim_MD_p012_2_res_to_glob",xlim,ids,1:3)

plot(xs, hcat(ys_cN,ys_work2),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3)))
plot(xs, hcat(ys_cN,ys_res2),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3)))
plot(xs, hcat(ys_cN,ys_w2glob),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3)))
plot(xs, hcat(ys_cN,ys_r2glob),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3)))



ys_cN_h16 = read_sums("test1_hh_16/sim_MD_p010_constN",xlim,ids,1:3)
ys_work4 = read_sums("test1_hh_16/sim_MD_p010_4_workloc",xlim,ids,1:3)
ys_res4 = read_sums("test1_hh_16/sim_MD_p010_4_resloc",xlim,ids,1:1)
ys_w4glob = read_sums("test1_hh_16/sim_MD_p010_4_work_to_glob",xlim,ids,1:1)
ys_r4glob = read_sums("test1_hh_16/sim_MD_p010_4_res_to_glob",xlim,ids,1:1)

plot(xs, hcat(ys_cN_h16,ys_work4),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3)))
plot(xs, hcat(ys_cN_h16,ys_res4),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 1)))
plot(xs, hcat(ys_cN_h16,ys_w4glob),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 1)))
plot(xs, hcat(ys_cN_h16,ys_r4glob),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 1)))


0


## MD p=0.15
xlim = floor(Int, 600/r_inc)
xs = collect(1:xlim) .* r_inc;
ys_base = read_sums("sim_MD_p015_base",xlim,ids,1:3)
ys_no_inc = read_sums("network_tests/sim_MD_p015_NO_INC",xlim,ids,1:3)
ys_r = read_sums("sim_MD_p015_rnet",xlim,ids,1:2)
ys_sw = read_sums("sim_MD_p015_sw",xlim,ids,1:3)
ys_ba = read_sums("network_tests/sim_MD_p015_ba",xlim,ids,1:3)
ys_er = read_sums("network_tests/sim_MD_p015_er",xlim,ids,1:2)
plot(xs, hcat(ys_er,ys_base,ys_ba),legend=nothing,
    linecolor=hcat(fill(:red, 1, 2), 
    fill(:green, 1, 3),
    fill(:black, 1, 3)))


## MD I=10,0,0, p=0.15
xlim = floor(Int, 600/r_inc)
xs = collect(1:xlim) .* r_inc;
ys_base = read_sums("sim_MD_p015_I10",xlim,ids,1:3)
ys_no_inc = read_sums("sim_MD_p015_I10_no_inc",xlim,ids,1:3)
plot(xs, hcat(ys_base,ys_no_inc),legend=nothing,
    linecolor=hcat(fill(:red, 1, 3), fill(:green, 1, 3)))



0


## high transmission (p = 0.4)
#xlim = floor(Int, 200/r_inc)
#xs = collect(1:xlim) .* r_inc
#ys_base = read_logs("sim_MD_p040_base",xlim,ids,1:3)
#ys_r = read_logs("sim_MD_p040_rnet",xlim,ids,1:2)
#ys_rrr = read_logs("sim_MD_p040_rrand",xlim,ids,1:2)
#plot(xs, hcat(ys_base,ys_r),legend=nothing)
#plot(xs, hcat(ys_r),legend=nothing)
#plot(xs, hcat(ys_r,ys_rrr),legend=nothing)

## very high transmission (p = 0.8)
#xlim = floor(Int, 100/r_inc)
#xs = collect(1:xlim) .* r_inc
#ys_base = read_logs("sim_MD_p080_base",xlim,ids,1:1)
#ys_r = read_logs("sim_MD_p080_rnet",xlim,ids,1:1)
#plot(xs, hcat(ys_base,ys_r),legend=nothing)
#plot(xs, hcat(ys_r),legend=nothing)






g = watts_strogatz(100, 8, 0.25)
g = barabasi_albert(1000, 4)
using CommunityDetection
x = community_detection_bethe(g,6)
[length(findall(a->a==i, x)) for i in 1:maximum(x)]

## try SBM for community detection


r,c,v = findnz(sprand(1000,1000,0.1))
x = sparse(r,c,trues(length(v)))

function myshuffle(x::SparseMatrixCSC{Bool, Int64})
    idxs = axes(x,1)
    r = Int64[]
    c = Int64[]
    for i in axes(x,2)
        n = nnz(view(x,:,i))
        k = sample(idxs, n, replace=false)
        append!(r,k)
        append!(c,fill(i,n))
    end
    return sparse(r,c,trues(nnz(x)))
end

@elapsed y = myshuffle(x)

function theirshuffle(x::SparseMatrixCSC{Bool, Int64})
    y = spzeros(Bool,size(x))
    for i in axes(x,2)
        y[:,i] = shuffle(x[:,i])
    end
    return y
end

@elapsed y = theirshuffle(x)



## access columns, not rows
@elapsed v = full_net[:,2]

## these are fast
@elapsed findnz(full_net)
@elapsed findall(full_net)

## access columns, not rows
@elapsed [findnz(full_net[:,i]) for i in 1:100]
@elapsed [findall(full_net[:,i]) for i in 1:100]

## column sums
t = sum(test_net,dims=1);
thh = sum(test_hh_net,dims=1);
f = sum(full_net,dims=1);

#[x=>sum(t.==x) for x in 0:12]
#[x=>sum(thh.==x) for x in 0:12]
#[x=>sum(f.==x) for x in 0:12]

#d = diag(full_net)









unit_ids = [2,3,4]
inputs = modelInputs(unit_ids; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI);
id = 3
in_chans = Dict(i => RemoteChannel(()->Channel{Vector{simEvent}}(1), i) for i in unit_ids)
out_chans = Dict(i => RemoteChannel(()->Channel{Vector{simEvent}}(1), i) for i in unit_ids)
report_chans = Dict(i => RemoteChannel(()->Channel{Any}(100), i) for i in unit_ids)
u = simUnit(in_chans[id], out_chans[id], report_chans[id])
u[:id] = id ## just for debugging?
init_sim_unit!(u, id, inputs) ## this fn adds domain-specific data and queues initial events

e = becomeContagious(3, 1000)
i = e.agentid

units = Dict()
for id in unit_ids
    u = simUnit(in_chans[id], out_chans[id], report_chans[id])
    u[:id] = id ## just for debugging?
    init_sim_unit!(u, id, inputs) ## this fn adds domain-specific data and queues initial events
    units[id] = u
end

unit = units[4];
for i in 1:50
    e, t = dequeue_pair!(unit[:q])	
    ## handle event; can change state, should generate and return future event(s)
    #println("unit ", unit[:id], " handling ", e)
    if ((e isa infectionEvent) || (e isa becomeContagious) || (e isa becomeRecovered))
        future_events = handle_event!(unit, e)
        sort_events!(unit, e, future_events)
    end
end

