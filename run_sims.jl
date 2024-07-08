using Distributed
addprocs(3; lazy=false)

include("gabm.jl")
include("sim_logic.jl")


## before first run:
println("precalc sets")
precalc_sets()
GC.gc()
##


unit_ids = workers()
n_units = length(unit_ids)
println("unit_ids = ",unit_ids)

function test(t; kwargs...)
    x = run(t, workers(); kwargs...)
    return x
end

t = 602
I0 = [2=>100,3=>100,4=>100]
pI = 0.1
t_closed = 50:400

outdir = "sim_close_wp"
nreps = 40

mkpath(outdir)

fname = "sim_p0"*string(round(Int,pI*100))*"_open"
for i in 1:nreps
    println("running ",fname," ",i,"/",nreps)
    serialize(outdir*"/"*fname*string(i)*".jlse", 
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=2.0*pI, p_inf_loc=pI,
    distr_fn_hh=:all, distr_fn_non_hh=:all, distr_params_hh=(), distr_params_non_hh=(), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(4,), distr_fn_loc_work=:exp, distr_params_loc_work=(2,) 
    ))
    GC.gc()
end

fname = "sim_p0"*string(round(Int,pI*100))*"_wp_sch_closed"
for i in 1:nreps
    println("running ",fname," ",i,"/",nreps)
    ## a random set of essential workers for each run
    ser_path("precalc/essential_workers.jlse", generate_ess_workers())
    serialize(outdir*"/"*fname*string(i)*".jlse",
    test(t; init_inf=I0, p_inf=pI, p_inf_hh=2.0*pI, p_inf_loc=pI, nonessential_wp_closed=t_closed, sch_closed=t_closed,
    distr_fn_hh=:all, distr_fn_non_hh=:all, distr_params_hh=(), distr_params_non_hh=(), 
    distr_fn_loc_res=:exp, distr_params_loc_res=(4,), distr_fn_loc_work=:exp, distr_params_loc_work=(2,) 
    ))
    GC.gc()
end

