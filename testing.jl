

include("gabm.jl")
include("sim_logic.jl")

##
## before first run:
##
## precalc_sets()


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
pI = 0.1
t_closed = 50:400

mkpath("test")

for i in 20:0
    serialize("test/sim_MD_p0"*string(round(Int,pI*100))*"_open"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=2.0*pI, p_inf_loc=pI,
    distr_fn_hh=:all, distr_fn_non_hh=:all, distr_params_hh=(), distr_params_non_hh=(), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(4,), distr_fn_loc_work=:exp, distr_params_loc_work=(2,) 
    ))
    GC.gc()
end


for i in 1:0
    serialize("test/sim_MD_p0"*string(round(Int,pI*100))*"_wp_closed"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=2.0*pI, p_inf_loc=pI, nonessential_wp_closed=t_closed,
    distr_fn_hh=:all, distr_fn_non_hh=:all, distr_params_hh=(), distr_params_non_hh=(), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(4,), distr_fn_loc_work=:exp, distr_params_loc_work=(2,) 
    ))
    GC.gc()
end

for i in 1:0
    ## a random set of essential workers for each run
    ser_path("precalc/essential_workers.jlse", generate_ess_workers())
    serialize("test/sim_MD_p0"*string(round(Int,pI*100))*"_wp_sch_closed"*string(i)*".jlse",
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=2.0*pI, p_inf_loc=pI, nonessential_wp_closed=t_closed, sch_closed=t_closed,
    distr_fn_hh=:all, distr_fn_non_hh=:all, distr_params_hh=(), distr_params_non_hh=(), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(4,), distr_fn_loc_work=:exp, distr_params_loc_work=(2,) 
    ))
    GC.gc()
end


#t = 10
@elapsed x = test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(8,), 
distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:zero, distr_params_loc_work=() 
);
for i in 6:0
    serialize("test/sim_MD_p0"*string(round(Int,pI*100))*"_base"*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
    distr_fn_hh=:const, distr_fn_non_hh=:const, distr_params_hh=(6,), distr_params_non_hh=(8,), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(2,), distr_fn_loc_work=:zero, distr_params_loc_work=() 
    ))
    GC.gc()
end


using Plots
using Shapefile

## read logs; arrange time series in matrix columns for Plots.plot()
#old_read_log(f::String,rep::Int) = deserialize(f*string(rep)*".jlse")[3]
## stack reps, sum workers
#old_read_sums(f,xlim,ids,reps) = hcat([sum(reduce(hcat, [[t[2] for t in old_read_log(f,r)[i][1:xlim]] for i in ids]),dims=2) for r in reps]...)

read_log(f::String,rep::Int) = deserialize(f*string(rep)*".jlse")
## time series sums:
read_sums(f,xlim,ids,k,reps) = reduce(hcat, sum([read_log(f,r)[i][k][1:xlim] for i in ids]) for r in reps)
## other sums:
nont_sums(f,ids,k,reps) = reduce(hcat, sum([read_log(f,r)[i][k] for i in ids]) for r in reps)

plotk(ns::Int,nr::Int) = [:legend=>nothing, :bg=>"#e1e4e1",
    :lc=>repeat(["#3b4b59" "#ed8008" "#736b1e" "#bf1b1b"][:,1:ns], inner=(1,nr)), 
    :ls=>repeat([:dashdot :solid :dash :dashdot][:,1:ns], inner=(1,nr))]
P(a,b)::Float64 = iszero(b) ? 0.0 : a/b

r_inc = 5
ids = 2:4
xlim = floor(Int, 600/r_inc)
xs = collect(1:xlim) .* r_inc;
test_time = Int(maximum(t_closed)/r_inc)

cbg_idxs = let k = dser_path("jlse/cbg_idxs.jlse"); Dict(String(v)=>Int(i) for (i,v) in k) end
pop_by_cbg = let d = dser_path("precalc/p_idxs_all_by_h_cbg.jlse"); Dict(String(k)=>length(v) for (k,v) in d) end
pop_by_cbgidx = Dict(v=>pop_by_cbg[k] for (k,v) in cbg_idxs)
pops_all = [pop_by_cbgidx[i] for i in 1:maximum(values(cbg_idxs))]


g_open = read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_open",xlim,ids,:geo_agents_assigned,1:20);
g_wp_sch = read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_wp_sch_closed",xlim,ids,:geo_agents_assigned,1:20);

total_open = read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_open",xlim,ids,:cumI,1:20);
total_wp_sch = read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_wp_sch_closed",xlim,ids,:cumI,1:20);
## time when open has the same cumI
test_time_o = argmin(abs.(vec(mean(total_open;dims=2)) .- mean(total_wp_sch[test_time,:])))



g_open_t = g_open[test_time_o,:]
g_open_p = [g_open_t[r] ./ pops_all for r in eachindex(g_open_t)]

g_wp_sch_t = g_wp_sch[test_time,:]
g_wp_sch_p = [g_wp_sch_t[r] ./ pops_all for r in eachindex(g_wp_sch_t)]


## baltimore
(shapes, cbg_codes) = let 
    t = Shapefile.Table("tl_2019_24_bg/tl_2019_24_bg.shp"); 
    mask = (t.COUNTYFP .== "510")
    (t.geometry[mask], t.GEOID[mask])
end;

## dc
(shapes, cbg_codes) = let 
    t = Shapefile.Table("tl_2019_11_bg/tl_2019_11_bg.shp")
    (t.geometry, t.GEOID)
end;


idxs = [get(cbg_idxs,k,-1) for k in cbg_codes]
#pops = [get(pop_by_cbg,k,0) for k in cbg_codes]


#g = mean(g_open_p); col_range = (0,0.4)
#g = sqrt.(var(g_open_p)); col_range = (0,0.1)
g = mean(g_wp_sch_p); col_range = (0,0.4)
#g = sqrt.(var(g_wp_sch_p)); col_range = (0,0.1)
z = reshape([get(g,i,0) for i in idxs], 1, :)


let k=:inferno; plot(shapes, clims=col_range, fill=palette(k), fill_z=z, linecolor=palette(k), line_z=z, size=(800,800)) end
#let k=:inferno; plot(shapes, fill=palette(k), fill_z=z, linecolor=palette(k), line_z=z, size=(800,800)) end





pidxs_by_cbg = dser_path("precalc/p_idxs_all_by_h_cbg.jlse");
n_b = let
    race_black = dser_path("precalc/race_black_alone.jlse")
    race_black_by_cbg = Dict(k => count(i->in(i,race_black), v) for (k,v) in pidxs_by_cbg)
    b_by_cbgidx = Dict(cbg_idxs[k]=>race_black_by_cbg[k] for k in keys(cbg_idxs))
    [b_by_cbgidx[i] for i in 1:maximum(values(cbg_idxs))] 
end
n_w =  let
    race_white = dser_path("precalc/white_non_hispanic.jlse")
    race_white_by_cbg = Dict(k => count(i->in(i,race_white), v) for (k,v) in pidxs_by_cbg)
    w_by_cbgidx = Dict(cbg_idxs[k]=>race_white_by_cbg[k] for k in keys(cbg_idxs))
    [w_by_cbgidx[i] for i in 1:maximum(values(cbg_idxs))] 
end
n_h = let
    hispanic_set = dser_path("precalc/hispanic.jlse")
    hispanic_by_cbg = Dict(k => count(i->in(i,hispanic_set), v) for (k,v) in pidxs_by_cbg)
    h_by_cbgidx = Dict(cbg_idxs[k]=>hispanic_by_cbg[k] for k in keys(cbg_idxs))
    [h_by_cbgidx[i] for i in 1:maximum(values(cbg_idxs))]
end


## baltimore
(shapes, cbg_codes) = let 
    t = Shapefile.Table("tl_2019_24_bg/tl_2019_24_bg.shp"); 
    mask = (t.COUNTYFP .== "510")
    (t.geometry[mask], t.GEOID[mask])
end;

## dc
(shapes, cbg_codes) = let 
    t = Shapefile.Table("tl_2019_11_bg/tl_2019_11_bg.shp")
    (t.geometry, t.GEOID)
end;


idxs = filter(i->i>0, [get(cbg_idxs,k,-1) for k in cbg_codes])

cI_b = read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_wp_sch_closed",xlim,ids,:geo_race_black_alone,1:20)[test_time,:];
cI_w = read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_wp_sch_closed",xlim,ids,:geo_white_non_hispanic,1:20)[test_time,:];
cI_h = read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_wp_sch_closed",xlim,ids,:geo_hispanic,1:20)[test_time,:];

pI_b = [sum(cI_b[r][idxs]) / sum(n_b[idxs]) for r in eachindex(cI_b)]
pI_w = [sum(cI_w[r][idxs]) / sum(n_w[idxs]) for r in eachindex(cI_w)]
pI_h = [sum(cI_h[r][idxs]) / sum(n_h[idxs]) for r in eachindex(cI_h)]

mean(pI_b), sqrt(var(pI_b))
mean(pI_w), sqrt(var(pI_w))
mean(pI_h), sqrt(var(pI_h))

24.06 / 19.0
(24.06 - 19.0) / 19.0
29.2 / 19.0
(29.2 - 19.0) / 19.0



cI_b = read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_open",xlim,ids,:geo_race_black_alone,1:20)[test_time_o,:];
cI_w = read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_open",xlim,ids,:geo_white_non_hispanic,1:20)[test_time_o,:];
cI_h = read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_open",xlim,ids,:geo_hispanic,1:20)[test_time_o,:];

pI_b = [sum(cI_b[r][idxs]) / sum(n_b[idxs]) for r in eachindex(cI_b)]
pI_w = [sum(cI_w[r][idxs]) / sum(n_w[idxs]) for r in eachindex(cI_w)]
pI_h = [sum(cI_h[r][idxs]) / sum(n_h[idxs]) for r in eachindex(cI_h)]

mean(pI_b), sqrt(var(pI_b))
mean(pI_w), sqrt(var(pI_w))
mean(pI_h), sqrt(var(pI_h))


25.8 / 22.3
(25.8 - 22.3) / 22.3
29.3 / 22.3
(29.3 - 22.3) / 22.3




cI_b = mean(read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_wp_sch_closed",xlim,ids,:geo_race_black_alone,1:20);dims=2)
cI_w = mean(read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_wp_sch_closed",xlim,ids,:geo_white_non_hispanic,1:20);dims=2)
cI_h = mean(read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_wp_sch_closed",xlim,ids,:geo_hispanic,1:20);dims=2)

cI_b_p = [sum(cI_b[i][idxs]) for i in eachindex(cI_b)] ./ sum(n_b[idxs])
cI_w_p = [sum(cI_w[i][idxs]) for i in eachindex(cI_w)] ./ sum(n_w[idxs])
cI_h_p = [sum(cI_h[i][idxs]) for i in eachindex(cI_h)] ./ sum(n_h[idxs])

plot(hcat((cI_b_p ./ cI_w_p), (cI_h_p ./ cI_w_p)))


cI_b = mean(read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_open",xlim,ids,:geo_race_black_alone,1:20);dims=2)
cI_w = mean(read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_open",xlim,ids,:geo_white_non_hispanic,1:20);dims=2)
cI_h = mean(read_sums("test/sim_MD_p0"*string(round(Int,pI*100))*"_open",xlim,ids,:geo_hispanic,1:20);dims=2)

cI_b_p = [sum(cI_b[i][idxs]) for i in eachindex(cI_b)] ./ sum(n_b[idxs])
cI_w_p = [sum(cI_w[i][idxs]) for i in eachindex(cI_w)] ./ sum(n_w[idxs])
cI_h_p = [sum(cI_h[i][idxs]) for i in eachindex(cI_h)] ./ sum(n_h[idxs])


plot(hcat((cI_b_p ./ cI_w_p), (cI_h_p ./ cI_w_p)))




## dc
t = Shapefile.Table("tl_2019_11_bg/tl_2019_11_bg.shp")
shapes = t.geometry
cbg_codes = t.GEOID

pops = [get(pop_by_cbg,k,0) for k in cbg_codes]
rblack = [get(race_black_by_cbg,k,0) for k in cbg_codes]
pI = reshape(P.(rblack,pops),1,:)
let z=pI; k=:inferno; plot(shapes, fill=palette(k), fill_z=z, linecolor=palette(k), line_z=z, size=(800,800)) end









## initial infections

p_by_cbg = dser_path("precalc/p_idxs_all_by_h_cbg.jlse")
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











## network tests
## MD p=0.15
xlim = floor(Int, 600/r_inc)
xs = collect(1:xlim) .* r_inc;
#ys_base = old_read_sums("bak/sim_MD_p015_base",xlim,ids,1:3)
#ys_no_inc = read_sums("network_tests/sim_MD_p015_NO_INC",xlim,ids,1:3)
#ys_r = old_read_sums("bak/sim_MD_p015_rnet",xlim,ids,1:3)
#ys_sw = old_read_sums("bak/sim_MD_p015_sw",xlim,ids,1:3)
#ys_ba = old_read_sums("bak/sim_MD_p015_ba",xlim,ids,1:3)
#ys_er = old_read_sums("bak/sim_MD_p015_er",xlim,ids,1:2)

xlim = 72
xs = collect(1:xlim) .* r_inc;
ys_base = old_read_sums("bak/sim_MD_base_netw",xlim,ids,1:3)
ys_ba = old_read_sums("bak/sim_MD_ba",xlim,ids,1:3)
ys_r = old_read_sums("bak/sim_MD_randonet",xlim,ids,1:2)
ys_rr = old_read_sums("bak/sim_MD_rrand",xlim,ids,1:2)

#3b4b59" "#ed8008" "#736b1e" "#bf1b1b
ns=3; nr=3;
plot(xs, [ys_ba ys_base  ys_r ys_rr[:,1]],
    label=["scale-free" "" "" "synth pop" "" "" "random network" "" ""],
    legend=:outerbottom,legendcolumns=3; 
    :bg=>"#e1e4e1",
    :lc=>repeat(["#3b4b59" "#ed8008" "#736b1e" "#bf1b1b"][:,1:ns], inner=(1,nr)), 
    :ls=>repeat([:dashdot :solid :dash :dashdot][:,1:ns], inner=(1,nr))
)

## MD I=10,0,0, p=0.15
xlim = floor(Int, 600/r_inc)
xs = collect(1:xlim) .* r_inc;
ys_base = read_sums("sim_MD_p015_I10",xlim,ids,1:3)
ys_no_inc = read_sums("sim_MD_p015_I10_no_inc",xlim,ids,1:3)
plot(xs, hcat(ys_base,ys_no_inc),legend=nothing,
    linecolor=hcat(fill(:red, 1, 3), fill(:green, 1, 3)))



0



## debugging event handler code



inputs = modelInputs(unit_ids; 
    init_inf=I0, p_inf=pI, p_inf_hh=2.0*pI, p_inf_loc=pI, nonessential_wp_closed=t_closed, sch_closed=t_closed,
    distr_fn_hh=:all, distr_fn_non_hh=:all, distr_params_hh=(), distr_params_non_hh=(), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(4,), distr_fn_loc_work=:exp, distr_params_loc_work=(2,) 
    #init_inf=I0, p_inf=pI, p_inf_hh=pI, p_inf_loc=pI,
    #flags=[:w_test],
    #distr_params_hh=(8,), distr_params_non_hh=(6,), distr_fn_loc_work=:const, distr_params_loc_work=(2,)
    );


in_chans = Dict(i => RemoteChannel(()->Channel{Vector{simEvent}}(1), i) for i in unit_ids)
out_chans = Dict(i => RemoteChannel(()->Channel{Vector{simEvent}}(1), i) for i in unit_ids)
report_chans = Dict(i => RemoteChannel(()->Channel{Any}(100), i) for i in unit_ids)


units = Dict()
for id in unit_ids
    u = SimUnit(id, in_chans[id], out_chans[id], report_chans[id])
    init_sim_unit!(u, inputs) ## this fn adds domain-specific data and queues initial events
    units[id] = u
end
#global_data = globalData(inputs)
inputs = nothing ; GC.gc()
#unit = units[3];
#q_event!(unit, syncEvent(unit.t_inc))
#t = UInt32(0)
for (k,unit) in units
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
end




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









using Clustering
assign_dict = Dict(i => UInt32[] for i in unit_ids)

let
    res_idxs_by_loc = dser_path("precalc/p_idxs_all_by_h_cbg.jlse")
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


