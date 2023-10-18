include("gabm.jl")
include("sim_logic.jl")

unit_ids = workers()
n_units = length(unit_ids)

function test(t; kwargs...)
    x = run(t, workers(); kwargs...)
    return x
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


for i in 6:6
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_base"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
    distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(8,), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:zero, distr_params_loc_work=() 
    ))
    GC.gc()
end




for i in 1:1
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_open"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
    distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(6,), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:exp, distr_params_loc_work=(2,) 
    ))
    GC.gc()
end

for i in 1:1
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_wp_closed"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI, nonessential_wp_closed=100:400,
    distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(6,), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:exp, distr_params_loc_work=(2,) 
    ))
    GC.gc()
end

for i in 1:1
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_wp_sch_closed"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI, nonessential_wp_closed=100:400, sch_closed=100:400,
    distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(6,), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:exp, distr_params_loc_work=(2,) 
    ))
    GC.gc()
end



@elapsed x = test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(8,), 
distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:zero, distr_params_loc_work=() 
)


using Plots

## read logs; arrange time series in matrix columns for Plots.plot()
old_read_log(f::String,rep::Int) = deserialize(f*string(rep)*".jlse")[3]
## stack reps, sum workers
old_read_sums(f,xlim,ids,reps) = hcat([sum(reduce(hcat, [[t[2] for t in old_read_log(f,r)[i][1:xlim]] for i in ids]),dims=2) for r in reps]...)

read_log(f::String,rep::Int) = deserialize(f*string(rep)*".jlse")
read_sums(f,xlim,ids,reps) = hcat([sum(reduce(hcat, [read_log(f,r)[i][:active][1:xlim] for i in ids]),dims=2) for r in reps]...)



r_inc = 5
ids = 2:4

xlim = floor(Int, 600/r_inc)
xs = collect(1:xlim) .* r_inc;


#ys_old = old_read_sums("sim_MD_p012_base",xlim,ids,1:1)
ys_base = read_sums("sim_MD_p012_base",xlim,ids,4:6)

plot(xs, ys_base, legend=nothing)
#plot(xs, [ys_old ys_base], legend=nothing, linecolor=[fill(:black, 1, 1) fill(:red, 1, 2)])






for i in 2:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_h152"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI, holiday_time=152:173,
    distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(8,), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:zero, distr_params_loc_work=() 
    ))
    GC.gc()
end

for i in 1:1
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_h348"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI, holiday_time=348:369,
    distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(8,), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:zero, distr_params_loc_work=() 
    ))
    GC.gc()
end


## initial infections

p_by_cbg = dser_path("p_idxs_all_by_h_cbg.jlse")
I0 = UInt32.(first(p_by_cbg["240230006001"], 100))

for i in 2:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_iloc"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
    distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(8,), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:zero, distr_params_loc_work=() 
    ))
    GC.gc()
end

for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_h202_iloc"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI, holiday_time=202:223,
    distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(8,), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:zero, distr_params_loc_work=() 
    ))
    GC.gc()
end



## testing local contacts

for i in 2:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_constN"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
    distr_params_hh=(8,)
    ))
    GC.gc()
end

for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2_workloc"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        flags=[:w_test],
        distr_params_hh=(8,), distr_params_non_hh=(6,), distr_fn_loc_work=:const, distr_params_loc_work=(2,)
        ))
    GC.gc()
end

for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2_resloc"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        flags=[:h_test],
        distr_params_hh=(6,), distr_fn_loc_res=:const, distr_params_loc_res=(2,)
        ))
    GC.gc()
end

for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2_res_to_glob"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        flags=[:h_test],
        distr_params_hh=(6,), distr_fn_nonloc=:const, distr_params_nonloc=(2,)
        ))
    GC.gc()
end

for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2_work_to_glob"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        flags=[:w_test],
        distr_params_hh=(8,), distr_params_non_hh=(6,), distr_fn_nonloc=:const, distr_params_nonloc=(2,)
        ))
    GC.gc()
end

for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2exp_workloc"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        flags=[:w_test],
        distr_params_hh=(8,), distr_params_non_hh=(6,), distr_fn_loc_work=:exp, distr_params_loc_work=(2,)
        ))
    GC.gc()
end

for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2exp_resloc"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        flags=[:h_test],
        distr_params_hh=(6,), distr_fn_loc_res=:exp, distr_params_loc_res=(2,)
        ))
    GC.gc()
end

for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2exp_res_to_glob"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        flags=[:h_test],
        distr_params_hh=(6,), distr_fn_nonloc=:exp, distr_params_nonloc=(2,)
        ))
    GC.gc()
end

for i in 1:3
    serialize("sim_MD_p0"*string(round(Int,pI*100))*"_2exp_work_to_glob"*string(i)*".jlse", 
        test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
        flags=[:w_test],
        distr_params_hh=(8,), distr_params_non_hh=(6,), distr_fn_nonloc=:exp, distr_params_nonloc=(2,)
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











## holiday tests
ys_h1 = read_sums("sim_MD_p012_h152",xlim,ids,1:3)
ys_h2 = read_sums("sim_MD_p012_h348",xlim,ids,1:1)
ys_iloc = read_sums("sim_MD_p012_iloc",xlim,ids,1:3)
ys_h1_iloc = read_sums("sim_MD_p012_h202_iloc",xlim,ids,1:3)
plot(xs, [ys_base ys_h1], legend=nothing, linecolor=[fill(:black, 1, 3) fill(:red, 1, 3)])
plot(xs, [ys_base ys_h2], legend=nothing, linecolor=[fill(:black, 1, 3) fill(:red, 1, 1)])
plot(xs, [ys_iloc ys_h1_iloc], legend=nothing, linecolor=[fill(:black, 1, 3) fill(:red, 1, 3)])

## location contacts
ys_cN = read_sums("sim_MD_p012_constN",xlim,ids,1:3)
ys_work2 = read_sums("sim_MD_p012_2_workloc",xlim,ids,1:3)
ys_res2 = read_sums("sim_MD_p012_2_resloc",xlim,ids,1:3)
ys_w2glob = read_sums("sim_MD_p012_2_work_to_glob",xlim,ids,1:3)
ys_r2glob = read_sums("sim_MD_p012_2_res_to_glob",xlim,ids,1:3)
ys_work2x = read_sums("sim_MD_p012_2exp_workloc",xlim,ids,1:3)
ys_res2x = read_sums("sim_MD_p012_2exp_resloc",xlim,ids,1:3)
ys_w2globx = read_sums("sim_MD_p012_2exp_work_to_glob",xlim,ids,1:3)
ys_r2globx = read_sums("sim_MD_p012_2exp_res_to_glob",xlim,ids,1:3)

plot(xs, hcat(ys_cN,ys_work2),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3)))
plot(xs, hcat(ys_cN,ys_res2),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3)))
plot(xs, hcat(ys_cN,ys_w2glob),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3)))
plot(xs, hcat(ys_cN,ys_r2glob),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3)))
plot(xs, hcat(ys_cN,ys_work2,ys_work2x),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3), fill(:blue, 1, 3)), linestyle=hcat(fill(:dash, 1, 3), fill(:solid, 1, 3), fill(:solid, 1, 3)))
plot(xs, hcat(ys_cN,ys_res2,ys_res2x),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3), fill(:blue, 1, 3)), linestyle=hcat(fill(:dash, 1, 3), fill(:solid, 1, 3), fill(:solid, 1, 3)))
plot(xs, hcat(ys_cN,ys_w2glob,ys_w2globx),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3), fill(:blue, 1, 3)), linestyle=hcat(fill(:dash, 1, 3), fill(:solid, 1, 3), fill(:solid, 1, 3)))
plot(xs, hcat(ys_cN,ys_r2glob,ys_r2globx),legend=nothing, linecolor=hcat(fill(:black, 1, 3), fill(:red, 1, 3), fill(:blue, 1, 3)), linestyle=hcat(fill(:dash, 1, 3), fill(:solid, 1, 3), fill(:solid, 1, 3)))

## network tests
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





## community detection?

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


##
## access columns, not rows (because Julia memory layout)
##
@elapsed v = full_net[:,2]

## these are fast (O1)
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







## debugging event handler code

inputs = modelInputs(unit_ids; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
    flags=[:w_test],
    distr_params_hh=(8,), distr_params_non_hh=(6,), distr_fn_loc_work=:const, distr_params_loc_work=(2,));
in_chans = Dict(i => RemoteChannel(()->Channel{Vector{simEvent}}(1), i) for i in unit_ids)
out_chans = Dict(i => RemoteChannel(()->Channel{Vector{simEvent}}(1), i) for i in unit_ids)
report_chans = Dict(i => RemoteChannel(()->Channel{Any}(100), i) for i in unit_ids)


units = Dict()
for id in unit_ids
    u = SimUnit(id, in_chans[id], out_chans[id], report_chans[id])
    init_sim_unit!(u, inputs) ## this fn adds domain-specific data and queues initial events
    units[id] = u
end

inputs = nothing ; GC.gc()

unit = units[4];
for i in 1:50
    e, t = dequeue_pair!(unit.q)	
    ## handle event; can change state, should generate and return future event(s)
    #println("unit ", unit[:id], " handling ", e)
    #if ((e isa infectionEvent) || (e isa becomeContagious) || (e isa becomeRecovered))
    if (!(e isa syncEvent))
        future_events = handle_event!(unit, e)
        sort_events!(unit, e, future_events)
    end
end






## graph partitioning by home cbg
## (this decreases # messages passed, but overall performance is the same b/c partitions not balanced)
## (balanced graph partitioning is NP-hard; but see approx algo Andreev 2004)


hh = dser_path("jlse/hh.jlse") ## households/residents (with hh locations)
gqs = dser_path("jlse/gqs.jlse") ## group-quarters/residents (with gq locations)
cbg_k = dser_path("jlse/cbg_idxs.jlse") ## location (cbg) keys used in person/hh keys
p_keys = dser_path("jlse/adj_mat_keys.jlse")
p_idxs = Dict(p_keys .=> eachindex(p_keys))

h_df_by_loc = groupby(DataFrame((k[2], v.people) for (k,v) in hh), "1")
hh_ppl_by_loc = Dict(loc["1"]=>reduce(vcat, h_df_by_loc[loc][!,"2"]) for loc in keys(h_df_by_loc))
hh_idxs_by_loc = Dict(k=>[p_idxs[i] for i in v] for (k,v) in hh_ppl_by_loc)

gq_df_by_loc = groupby(DataFrame((k[2], v.residents) for (k,v) in gqs), "1")
gq_ppl_by_loc = Dict(loc["1"]=>reduce(vcat, gq_df_by_loc[loc][!,"2"]) for loc in keys(gq_df_by_loc))
gq_idxs_by_loc = Dict(k=>[p_idxs[i] for i in v] for (k,v) in gq_ppl_by_loc)

## residential locations include households and non-inst GQs:
res_idxs_by_loc = vecmerge(hh_idxs_by_loc, gq_idxs_by_loc)

sum(length.(values(res_idxs_by_loc)))

work_dummies = dser_path("jlse/work_dummies.jlse")
dummy_idxs = [p_idxs[k[1:3]] for k in work_dummies]
res_idxs_by_loc[0] = dummy_idxs

sum(length.(values(res_idxs_by_loc)))

cbg_k[0] = "outside"
res_idxs_by_loc = Dict(cbg_k[k]=>v for (k,v) in res_idxs_by_loc)

ser_path("p_idxs_all_by_h_cbg.jlse", res_idxs_by_loc)




using Clustering
assign_dict = Dict(i => UInt32[] for i in unit_ids)

let
    res_idxs_by_loc = dser_path("p_idxs_all_by_h_cbg.jlse")
    d = read_df("../sim-netw/processed/work_od_matrix_no_inc.csv";types=Dict("h_cbg"=>String15))
    home_labels = d[:,1]
    M = Matrix{Float64}(d[!,2:end])
    M = Matrix(M') ## home cbg as columns
    result = kmeans(M, n_units; maxiter=500, tol=1.0e-7)

    cluster_idx = Dict(home_labels .=> result.assignments)
    cluster_idx["outside"] = argmin(result.counts) ## put these wherever for now
    for (k,v) in res_idxs_by_loc
        append!(assign_dict[unit_ids[cluster_idx[k]]], v)
    end
end

assign_dict
[(k,length(v)) for (k,v) in assign_dict]



using Distributions
using StatsPlots

p2 = Poisson(2)
mean(p2)
median(p2)
std(p2)
boxplot(rand(p2,1000000))

e2 = Exponential(2.0)
mean(e2)
median(e2)
std(e2)
violin(rand(e2,1000000))
boxplot(rand(e2,100000))

g2 = Geometric(1/(1+2.0))
mean(g2)
median(g2)
std(g2)
boxplot([rand(e2,100000) rand(g2,100000)])



## make dicts callable
#(d::Dict)(x) = d[x]
#p_lookup = let t = dser_path("jlse/adj_mat_keys.jlse"); Dict(t .=> UInt32.(eachindex(t))) end
#ppl_by_hh = let t = dser_path("jlse/hh_ppl.jlse"); Dict(hkey => map(p_lookup, pkeys) for (hkey,pkeys) in t) end
#hh_lookup = Dict(v=>k[2:3] for (k,v) in p_lookup)









####
function OLD_handle_event!(u::SimUnit, e::becomeContagious)::Vector{simEvent}
    i = e.agentid
    if resists_infection(i, e, u)
        return simEvent[] ## agent not susceptible
    else
        u[:cum_I] += 1 ## update cumulative infection count
        push!(u[:I_set], i) ## append to current infected set
        duration = rand(u[:t_recovery])
        recov_event = becomeRecovered(e.t + duration, i)

        ## note, findall() on a sparsearray is not really a "find", it's an O(1) lookup
        neigh_hh::Vector{UInt32} = findall(u[:netw_hh][:,i]) 
        #neigh_non_hh::Vector{UInt32} = findall(u[:netw_non_hh][:,i]) 
        neigh_wp::Vector{UInt32} = findall(u[:netw_wp][:,i]) 
        neigh_sch::Vector{UInt32} = findall(u[:netw_sch][:,i]) 
        neigh_gq::Vector{UInt32} = findall(u[:netw_gq][:,i]) 
        neigh_non_hh = [neigh_wp; neigh_sch; neigh_gq]

        ## possible ephemeral/location-based contacts
        ## (matrix columns are locations)
        neigh_loc_res::Vector{UInt32} = haskey(u[:loc_lookup_res], i) ? findall(u[:loc_matrix_res][:,u[:loc_lookup_res][i]]) : UInt32[]
        neigh_loc_work::Vector{UInt32} = haskey(u[:loc_lookup_work], i) ? findall(u[:loc_matrix_work][:,u[:loc_lookup_work][i]]) : UInt32[]

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

        contacts_non_hh::Vector{UInt32} = u[:distr_fn_non_hh](neigh_non_hh, u[:distr_params_non_hh])
        contacts_hh::Vector{UInt32} = u[:distr_fn_hh](neigh_hh, u[:distr_params_hh])
        contacts_loc_res::Vector{UInt32} = u[:distr_fn_loc_res](neigh_loc_res, distr_params_loc_res)
        contacts_loc_work::Vector{UInt32} = u[:distr_fn_loc_work](neigh_loc_work, distr_params_loc_work)
        ## axes 1: anyone in the population is a potential nonlocal contact
        contacts_nonloc::Vector{UInt32} = u[:distr_fn_nonloc](axes(u[:netw_hh],1), distr_params_nonloc) 
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



## for now, just figure out how many people you would have infected, and make them random people instead
## (keep infection rate unchanged to test the effect of network disruption)
function OLD_handle_event!(u::SimUnit, e::becomeHolidayContagious)::Vector{simEvent}
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

        everyone = axes(u[:netw_non_hh],1)

        ## doing this because distr fn might not really be a distr fn
        contacts_non_hh::Vector{UInt32} = u[:distr_fn_non_hh](neigh_non_hh, u[:distr_params_non_hh])
        contacts_hh::Vector{UInt32} = u[:distr_fn_hh](neigh_hh, u[:distr_params_hh])
        contacts_loc_res::Vector{UInt32} = u[:distr_fn_loc_res](neigh_loc_res, u[:distr_params_loc_res])
        contacts_loc_work::Vector{UInt32} = u[:distr_fn_loc_work](neigh_loc_work, u[:distr_params_loc_work])
        contacts_nonloc::Vector{UInt32} = u[:distr_fn_nonloc](everyone, u[:distr_params_nonloc]) 

        ## doing this in case p_inf's are different
        contacts_non_hh = rsamp(everyone, length(contacts_non_hh))
        contacts_hh = rsamp(everyone, length(contacts_hh))
        contacts_loc_res = rsamp(everyone, length(contacts_loc_res))
        contacts_loc_work = rsamp(everyone, length(contacts_loc_work))
        contacts_nonloc = rsamp(everyone, length(contacts_nonloc))

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


function OLD_modelInputs(unit_ids::Vector{Int64}; kwargs...)

    ## these funcs only evaluated when needed
    read_dummies() = Set{UInt32}(keys(dser_path("jlse/adj_dummy_keys.jlse")))
    read_outw() = Set{UInt32}(keys(dser_path("jlse/adj_out_workers.jlse")))
    read_netw_hh() = SparseMatrixCSC{Bool,UInt32}(Symmetric(dser_path("jlse/adj_mat_hh.jlse")))
    #read_netw_non_hh() = SparseMatrixCSC{Bool,UInt32}(Symmetric(dser_path("jlse/adj_mat_non_hh.jlse")))
    read_netw_wp() = SparseMatrixCSC{Bool,UInt32}(Symmetric(dser_path("jlse/adj_mat_wp.jlse")))
    read_netw_sch() = SparseMatrixCSC{Bool,UInt32}(Symmetric(dser_path("jlse/adj_mat_sch.jlse")))
    read_netw_gq() = SparseMatrixCSC{Bool,UInt32}(Symmetric(dser_path("jlse/adj_mat_gq.jlse")))
    read_loc_matrix_res() = SparseMatrixCSC{Bool,UInt32}(dser_path("jlse/res_loc_contact_mat.jlse"))
    read_loc_matrix_work() = SparseMatrixCSC{Bool,UInt32}(dser_path("jlse/work_loc_contact_mat.jlse"))
    read_loc_lookup_res() = Dict{UInt32,UInt32}(dser_path("jlse/res_loc_lookup.jlse"))
    read_loc_lookup_work() = Dict{UInt32,UInt32}(dser_path("jlse/work_loc_lookup.jlse"))
    #read_netw_holiday() = SparseMatrixCSC{Bool,UInt32}(Symmetric(dser_path("jlse/holiday_adj_mat.jlse")))
    #read_loc_matrix_holiday() = SparseMatrixCSC{Bool,UInt32}(dser_path("jlse/holiday_loc_contact_mat.jlse"))
    #read_loc_lookup_holiday() = Dict{UInt32,UInt32}(dser_path("jlse/holiday_loc_lookup.jlse"))


    ## awkward syntax but avoids costly evaluation when kwargs has the key
    netw_hh::SparseMatrixCSC{Bool, UInt32} = get(read_netw_hh, kwargs, :netw_hh)
    #netw_non_hh::SparseMatrixCSC{Bool, UInt32} = get(read_netw_non_hh, kwargs, :netw_non_hh)
    netw_wp::SparseMatrixCSC{Bool, UInt32} = get(read_netw_wp, kwargs, :netw_wp)

    ## dummies appear in work networks, but live outside the synth pop (no household or demographic info generated)
    ## local sim unit is responsible for determining if they got infected at home
    dummies::Set{UInt32} = get(read_dummies, kwargs, :dummies)
    ## outside workers have households, but no workplace network
    ## local sim unit is responsible for determining if they got infected at work
    outw::Set{UInt32} = get(read_outw, kwargs, :out_workers)

    mu_hh_cnx::Float64 = get(()->calc_mean_hh(netw_hh,dummies), kwargs, :mean_hh_connections)
    #mu_work_cnx::Float64 = get(()->calc_mean_wp(netw_non_hh,outw), kwargs, :mean_wp_connections)
    mu_work_cnx::Float64 = get(()->calc_mean_wp(netw_wp,outw), kwargs, :mean_wp_connections)

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
        #netw_non_hh = netw_non_hh,
        netw_wp = netw_wp,
        netw_sch = get(read_netw_sch, kwargs, :netw_sch),
        netw_gq = get(read_netw_gq, kwargs, :netw_gq),
        loc_matrix_res = get(read_loc_matrix_res, kwargs, :loc_matrix_res),
        loc_matrix_work = get(read_loc_matrix_work, kwargs, :loc_matrix_work),
        loc_lookup_res = get(read_loc_lookup_res, kwargs, :loc_lookup_res),
        loc_lookup_work = get(read_loc_lookup_work, kwargs, :loc_lookup_work),
        #netw_holiday = get(read_netw_holiday, kwargs, :netw_holiday),
        #loc_matrix_holiday = get(read_loc_matrix_holiday, kwargs, :loc_matrix_holiday),
        #loc_lookup_holiday = get(read_loc_lookup_holiday, kwargs, :loc_lookup_holiday),
        #netw = get(()->calc_full_net(netw_hh,netw_non_hh), kwargs, :full_net),
        t_inc = get(kwargs, :t_inc, 5),
        t_recovery = get(kwargs, :t_recovery, 8:12),
        #holiday_time = get(kwargs, :holiday_time, 0:0),
        init_inf = get(kwargs, :init_inf, [first(unit_ids) => 10]),
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

