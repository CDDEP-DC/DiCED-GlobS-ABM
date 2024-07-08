include("fileutils.jl")
include("utils.jl")

using Plots
using Shapefile

read_log(f::String,rep::Int) = deserialize(f*string(rep)*".jlse")
## time series sums:
read_sums(f,xlim,ids,k,reps) = reduce(hcat, sum([read_log(f,r)[i][k][1:xlim] for i in ids]) for r in reps)
## other sums:
nont_sums(f,ids,k,reps) = reduce(hcat, sum([read_log(f,r)[i][k] for i in ids]) for r in reps)

plotk(ns::Int,nr::Int) = [:legend=>nothing, :bg=>"#e1e4e1",
    :lc=>repeat(["#3b4b59" "#ed8008" "#736b1e" "#bf1b1b"][:,1:ns], inner=(1,nr)), 
    :ls=>repeat([:dashdot :solid :dash :dashdot][:,1:ns], inner=(1,nr))]
P(a,b)::Float64 = iszero(b) ? 0.0 : a/b

function draw_geo(vals, map_idxs, shapes, name)
    z = reshape([get(vals,i,0) for i in map_idxs], 1, :)
    k=:inferno
    col_range = (0,1.75)
    plot(shapes, clims=col_range, fill=palette(k), fill_z=z, linecolor=palette(k), line_z=z, size=(800,800), dpi=300) 
    savefig(name*".png")
end

function plot_rr(rr_open,rr_closed,xs,ymax,name)
    plot(xs,[mean(rr_open;dims=2) mean(rr_closed;dims=2)], 
    linecolor=[:black :black], linestyle=[:dash :solid],
    ylim=(0.9,ymax),
    label=["open" "closed"],ylabel="relative risk",xlabel="day",
    dpi=300)

    plot!(xs,[percentile.(eachrow(rr_open),10) percentile.(eachrow(rr_closed),10)], 
    fillrange=[percentile.(eachrow(rr_open),90) percentile.(eachrow(rr_closed),90)],
    alpha=0.3,linewidth=0,fillcolor=[:gray55 :gray33],
    label="")

    savefig(name*".png")
end


cbg_idxs = let k = dser_path("jlse/cbg_idxs.jlse"); Dict(String(v)=>Int(i) for (i,v) in k) end
pop_by_cbg = let d = dser_path("precalc/p_idxs_all_by_h_cbg.jlse"); Dict(String(k)=>length(v) for (k,v) in d) end
pop_by_cbgidx = Dict(v=>pop_by_cbg[k] for (k,v) in cbg_idxs)
pops_all = [pop_by_cbgidx[i] for i in 1:maximum(values(cbg_idxs))]

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


r_inc = 5
ids = 2:4
xlim = floor(Int, 600/r_inc)
test_time = 30
pI = 0.1
outdir = "sim_close_wp"
nreps = 40
f_open = "sim_p0"*string(round(Int,pI*100))*"_open"
f_closed = "sim_p0"*string(round(Int,pI*100))*"_wp_sch_closed"

g_open_t = read_sums(outdir*"/"*f_open,xlim,ids,:geo_agents_assigned,1:nreps)[test_time,:];
g_wp_sch_t = read_sums(outdir*"/"*f_closed,xlim,ids,:geo_agents_assigned,1:nreps)[test_time,:];
g_open_p = [g_open_t[r] ./ pops_all for r in eachindex(g_open_t)]
g_wp_sch_p = [g_wp_sch_t[r] ./ pops_all for r in eachindex(g_wp_sch_t)]

total_p_open = [sum(g_open_t[r]) / sum(pops_all) for r in eachindex(g_open_t)]
total_p_wp_sch = [sum(g_wp_sch_t[r]) / sum(pops_all) for r in eachindex(g_wp_sch_t)]


dGeo = tryJSON("geo.json")
map_counties = get(dGeo, "draw_map", [])
stfps = unique([x[1:2] for x in map_counties])
tables = [Shapefile.Table("geo/tl_2019_"*x*"_bg.shp") for x in stfps];
shapes = reduce(vcat, [t.geometry for t in tables])
cbg_codes = reduce(vcat, [t.GEOID for t in tables])
cbg_county = [x[1:5] for x in cbg_codes]
mask = [x in map_counties for x in cbg_county]
idxs = [get(cbg_idxs,k,-1) for k in cbg_codes]
gname = join(map_counties,"_")
draw_geo(mean(g_open_p ./ total_p_open), idxs[mask], shapes[mask], "map_"*gname*"_open")
draw_geo(mean(g_wp_sch_p ./ total_p_wp_sch), idxs[mask], shapes[mask], "map_"*gname*"_closed")


## whole pop
cbg_codes = collect(keys(cbg_idxs))
idxs = filter(i->i>0, [get(cbg_idxs,k,-1) for k in cbg_codes])
gname = join(unique([x[1:2] for x in get(dGeo, "geos", [])]),"_")

f(v)=sum(v[idxs])
cI_b = f.(read_sums(outdir*"/"*f_closed,xlim,ids,:geo_race_black_alone,1:nreps));
cI_w = f.(read_sums(outdir*"/"*f_closed,xlim,ids,:geo_white_non_hispanic,1:nreps));
cI_h = f.(read_sums(outdir*"/"*f_closed,xlim,ids,:geo_hispanic,1:nreps));

oI_b = f.(read_sums(outdir*"/"*f_open,xlim,ids,:geo_race_black_alone,1:nreps));
oI_w = f.(read_sums(outdir*"/"*f_open,xlim,ids,:geo_white_non_hispanic,1:nreps));
oI_h = f.(read_sums(outdir*"/"*f_open,xlim,ids,:geo_hispanic,1:nreps));

cI_b_p = cI_b ./ sum(n_b[idxs])
cI_w_p = cI_w ./ sum(n_w[idxs])
cI_h_p = cI_h ./ sum(n_h[idxs])

oI_b_p = oI_b ./ sum(n_b[idxs])
oI_w_p = oI_w ./ sum(n_w[idxs])
oI_h_p = oI_h ./ sum(n_h[idxs])

rr_c_h = cI_h_p ./ cI_w_p
rr_c_b = cI_b_p ./ cI_w_p

rr_o_h = oI_h_p ./ oI_w_p
rr_o_b = oI_b_p ./ oI_w_p

xs = 5 .* collect(10:120)
plot_rr(rr_o_b[10:end,:], rr_c_b[10:end,:], xs, 1.55, "rr_"*gname*"_black")
plot_rr(rr_o_h[10:end,:], rr_c_h[10:end,:], xs, 2.05, "rr_"*gname*"_hispanic")





