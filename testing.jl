include("gabm.jl")
include("sim_logic.jl")

function test(t; kwargs...)
    x = run(t, workers(); kwargs...)
    return Dict(k => (v[:cum_I],v[:curr_I],v[:curr_R],v[:q_len],v[:glob_len]) for (k,v) in x),
            Dict(k => v[:n_agents] for (k,v) in x),
            Dict(k => v[:log] for (k,v) in x)
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
pI = 0.15

for i in 1:3
    serialize("sim_MD_p015_base"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI))
    GC.gc()
end

for i in 1:3
    randonet = colshuffle(sparse(Symmetric(dser_path("jlse/adj_mat.jlse"))) .| sparse(Symmetric(dser_path("jlse/hh_adj_mat.jlse"))));
    GC.gc()
    serialize("sim_MD_p015_rnet"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI, full_net=randonet))
    GC.gc()
end


#netw = sparse(Symmetric(dser_path("jlse/adj_mat.jlse"))) .| sparse(Symmetric(dser_path("jlse/hh_adj_mat.jlse")));
#dum_idxs = Set(keys(dser_path("jlse/adj_dummy_keys.jlse")));
#outw_idxs = Set(keys(dser_path("jlse/adj_out_workers.jlse")));
#mask = setdiff(axes(netw,2), dum_idxs, outw_idxs)
#mu = mean(sum(netw[:,mask], dims=1))
#nn = size(netw,1)
#netw, dum_idxs, outw_idxs, mask = nothing, nothing, nothing, nothing
#GC.gc()
#mu = 8.127311432232107 ## LA
#nn = 18299561 ## LA
mu = 8.050331099788291 ## MD
nn = 9440777 ## MD

for i in 1:3
    really_rand = randnet(nn, mu, 1.0)
    GC.gc()
    serialize("sim_MD_p015_rrand"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI, full_net=really_rand, dummies = Set(Int64[]), out_workers = Set(Int64[])))
    GC.gc()
end

using Graphs
nn = 9440777 ## MD

for i in 1:3
    sw_net = Bool.(Graphs.LinAlg.adjacency_matrix( watts_strogatz(nn,8,0.25) ));
    GC.gc()
    serialize("sim_MD_p015_sw"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI, full_net=sw_net, dummies = Set(Int64[]), out_workers = Set(Int64[])))
    GC.gc()
end

for i in 1:3
    ba_net = Bool.(Graphs.LinAlg.adjacency_matrix( barabasi_albert(nn, 4) ));
    GC.gc()
    serialize("sim_MD_p015_ba"*string(i)*".jlse", test(t; init_inf=I0, p_inf=pI, full_net=ba_net, dummies = Set(Int64[]), out_workers = Set(Int64[])))
    GC.gc()
end



using Plots

read_log(f::String,rep::Int) = deserialize(f*string(rep)*".jlse")[3]
read_logs(f,xlim,ids,reps) = hcat([[t[2] for t in read_log(f,rep)[i][1:xlim]] for i in ids for rep in reps]...)

r_inc = 5
ids = workers()
xlim = floor(Int,t/r_inc)
xs = collect(1:xlim) .* r_inc

## MD p 0.15
ys_base = read_logs("sim_MD_p015_base",xlim,ids,1:3)
ys_r = read_logs("sim_MD_p015_rnet",xlim,ids,1:3)
ys_rrr = read_logs("sim_MD_p015_rrand",xlim,ids,1:3)
ys_sw = read_logs("sim_MD_p015_sw",xlim,ids,1:3)
ys_ba = read_logs("sim_MD_p015_ba",xlim,ids,1:3)
plot(xs, hcat(ys_base,ys_r,ys_rrr),legend=nothing)
plot(xs, hcat(ys_base,ys_r,ys_rrr,ys_sw,ys_ba),legend=nothing)


xlim = 72
xs2 = collect(1:xlim) .* r_inc

## MD
ys_base2 = read_logs("sim_MD_base_netw",xlim,ids,1:3)
ys_r2 = read_logs("sim_MD_randonet",xlim,ids,1:2)
ys_rrr2 = read_logs("sim_MD_rrand",xlim,ids,1:2)
ys_sw2 = read_logs("sim_MD_sw",xlim,ids,1:3)
ys_ba2 = read_logs("sim_MD_ba",xlim,ids,1:3)
plot(xs2, hcat(ys_base2,ys_r2,ys_rrr2),legend=nothing)
plot(xs2, hcat(ys_base2,ys_r2,ys_rrr2,ys_sw2,ys_ba2[:,4:6]),legend=nothing)



## LA
ys_base = read_logs("sim_364_base_netw",xlim,ids,1:3)
ys_r = read_logs("sim_364_randonet",xlim,ids,1:3)
ys_rrr = read_logs("sim_364_rrand",xlim,ids,1:2)
plot(xs2, hcat(ys_base,ys_r,ys_rrr),legend=nothing)



g = watts_strogatz(100, 8, 0.25)
g = barabasi_albert(1000, 4)
using CommunityDetection
x = community_detection_bethe(g,6)
[length(findall(a->a==i, x)) for i in 1:maximum(x)]




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

