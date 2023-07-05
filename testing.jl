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
t = 364
I0 = [2=>100,3=>100,4=>100]


end_states, pop_sizes, logs = test(t; init_inf=I0)
#serialize("sim_MD_base_netw1.jlse", (end_states, pop_sizes, logs))
#serialize("sim_MD_base_netw2.jlse", (end_states, pop_sizes, logs))
#serialize("sim_364_base_netw1.jlse", (end_states, pop_sizes, logs))
#serialize("sim_364_base_netw2.jlse", (end_states, pop_sizes, logs))
#serialize("sim_364_base_netw3.jlse", (end_states, pop_sizes, logs))
GC.gc()


netw = sparse(Symmetric(dser_path("jlse/adj_mat.jlse"))) .| sparse(Symmetric(dser_path("jlse/hh_adj_mat.jlse")));
GC.gc()
randonet = colshuffle(netw)
netw = nothing
GC.gc()

end_states, pop_sizes, logs = test(t; init_inf=I0, full_net=randonet)

#serialize("sim_MD_randonet1.jlse", (end_states, pop_sizes, logs))
#serialize("sim_MD_randonet2.jlse", (end_states, pop_sizes, logs))
#serialize("sim_364_randonet1.jlse", (end_states, pop_sizes, logs))
#serialize("sim_364_randonet2.jlse", (end_states, pop_sizes, logs))
#serialize("sim_364_randonet3.jlse", (end_states, pop_sizes, logs))
GC.gc()



#netw = sparse(Symmetric(dser_path("jlse/adj_mat.jlse"))) .| sparse(Symmetric(dser_path("jlse/hh_adj_mat.jlse")));
#dum_idxs = Set(keys(dser_path("jlse/adj_dummy_keys.jlse")));
#outw_idxs = Set(keys(dser_path("jlse/adj_out_workers.jlse")));
#mask = setdiff(axes(netw,2), dum_idxs, outw_idxs)
#mu = mean(sum(netw[:,mask], dims=1))
#nv = size(netw,1)
#netw, dum_idxs, outw_idxs, mask = nothing, nothing, nothing, nothing
#GC.gc()
#mu = 8.127311432232107 ## LA
#nv = 18299561 ## LA

mu = 8.050331099788291 ## MD
nv = 9440777 ## MD
really_rand = randnet(nv, mu, 1.0)
GC.gc()
end_states, pop_sizes, logs = test(t; init_inf=I0, full_net=really_rand, dummies = Set(Int64[]), out_workers = Set(Int64[]))

#serialize("sim_MD_rrand1.jlse", (end_states, pop_sizes, logs))
#serialize("sim_MD_rrand2.jlse", (end_states, pop_sizes, logs))
#serialize("sim_364_rrand1.jlse", (end_states, pop_sizes, logs))
#serialize("sim_364_rrand2.jlse", (end_states, pop_sizes, logs))
#serialize("sim_364_rrand3.jlse", (end_states, pop_sizes, logs))



using Plots

end_states, pop_sizes, logs = deserialize("sim_364_base_netw1.jlse");
xs = [Int64(t[1]) for t in logs[2]];
ys_base1 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])
end_states, pop_sizes, logs = deserialize("sim_364_base_netw2.jlse");
ys_base2 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])
end_states, pop_sizes, logs = deserialize("sim_364_base_netw3.jlse");
ys_base3 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])

end_states, pop_sizes, logs = deserialize("sim_364_randonet1.jlse");
ys_r1 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])
end_states, pop_sizes, logs = deserialize("sim_364_randonet2.jlse");
ys_r2 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])
end_states, pop_sizes, logs = deserialize("sim_364_randonet3.jlse");
ys_r3 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])

end_states, pop_sizes, logs = deserialize("sim_364_rrand1.jlse");
ys_rrr1 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])
end_states, pop_sizes, logs = deserialize("sim_364_rrand2.jlse");
ys_rrr2 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])
end_states, pop_sizes, logs = deserialize("sim_364_rrand3.jlse");
ys_rrr3 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])

plot(xs, hcat(ys_base2,ys_r2,ys_rrr2))



using Plots

end_states, pop_sizes, logs = deserialize("sim_MD_base_netw1.jlse");
xs = [Int64(t[1]) for t in logs[2]];
ys_base1 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])
end_states, pop_sizes, logs = deserialize("sim_MD_base_netw2.jlse");
ys_base2 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])

end_states, pop_sizes, logs = deserialize("sim_MD_randonet1.jlse");
ys_r1 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])
end_states, pop_sizes, logs = deserialize("sim_MD_randonet2.jlse");
ys_r2 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])

end_states, pop_sizes, logs = deserialize("sim_MD_rrand1.jlse");
ys_rrr1 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])
end_states, pop_sizes, logs = deserialize("sim_MD_rrand2.jlse");
ys_rrr2 = reduce(hcat, [[t[2] for t in logs[i]] for i in 2:4])

plot(xs, hcat(ys_base1,ys_r1,ys_r2,ys_rrr1,ys_rrr2))






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

