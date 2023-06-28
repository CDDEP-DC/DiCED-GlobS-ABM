include("gabm.jl")
include("sim_logic.jl")

full_net = sparse(Symmetric(dser_path("jlse/adj_mat.jlse"))) .| sparse(Symmetric(dser_path("jlse/hh_adj_mat.jlse")));
GC.gc()

function test(t::Int)
    x = run(t, workers(); network=full_net, p_inf=0.2, init_inf=[2=>10]);
    return Dict(k => (v[:cum_I],v[:curr_I],v[:curr_R],v[:q_len],v[:glob_len]) for (k,v) in x)
end

test(300)


#netw = dser_path("jlse/netw.jlse")
#hh_net = dser_path("jlse/hh_ppl.jlse")
#owork = dser_path("jlse/outside_workers.jlse")
#cbgs = sort(unique(k[2] for k in keys(hh_net)))

## access columns, not rows
@elapsed v = full_net[:,2]

@elapsed findnz(full_net)
@elapsed findall(full_net)

## access columns, not rows
@elapsed [findnz(full_net[:,i]) for i in 1:100]
@elapsed [findall(full_net[:,i]) for i in 1:100]

t = sum(test_net,dims=1);
thh = sum(test_hh_net,dims=1);
f = sum(full_net,dims=1);

#[x=>sum(t.==x) for x in 0:12]
#[x=>sum(thh.==x) for x in 0:12]
#[x=>sum(f.==x) for x in 0:12]

#d = diag(full_net)

