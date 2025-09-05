
include("sim_undiced.jl")

## pre-generated ## precalc_params()
param_sets = dser_path("param_sets.jlse");

## keep in memory to avoid reloading each run
netw_hh = SparseMatrixCSC{Bool, UInt32}(Symmetric(dser_path("jlse/adj_mat_hh.jlse")));
agents_assigned = Set(UInt32.(collect(1:size(netw_hh,2))))
netw_wp = SparseMatrixCSC{Bool, UInt32}(Symmetric(dser_path("jlse/adj_mat_wp.jlse")));
dummies = collect(UInt32, keys(dser_path("jlse/adj_dummy_keys.jlse")));
out_workers = collect(UInt32, keys(dser_path("jlse/adj_out_workers.jlse")));
netw_sch = SparseMatrixCSC{Bool, UInt32}(Symmetric(dser_path("jlse/adj_mat_sch.jlse")));
netw_gq = SparseMatrixCSC{Bool, UInt32}(Symmetric(dser_path("jlse/adj_mat_gq.jlse")));
loc_matrix_res = SparseMatrixCSC{Bool,UInt32}(dser_path("jlse/res_loc_contact_mat.jlse"));
loc_matrix_work = SparseMatrixCSC{Bool,UInt32}(dser_path("jlse/work_loc_contact_mat.jlse"));
loc_lookup_res = Dict{UInt32,UInt32}(dser_path("jlse/res_loc_lookup.jlse"));
loc_lookup_work = Dict{UInt32,UInt32}(dser_path("jlse/work_loc_lookup.jlse"));
geo_lookup = Dict{UInt32,UInt16}(dser_path("precalc/cbg_idx_lookup.jlse"));
k12_workers = Set{UInt32}(dser_path("precalc/k12_workers.jlse"));
age_65o = Set{UInt32}(dser_path("precalc/age_65o.jlse"));
mean_hh_connections = Float64(calc_mean_hh(netw_hh,dummies))
mean_wp_connections = Float64(calc_mean_wp(netw_wp,out_workers))

refs = Dict(
        :netw_hh => netw_hh,
        :agents_assigned => agents_assigned,
        :netw_wp => netw_wp,
        :dummies => dummies,
        :out_workers => out_workers,
        :netw_sch => netw_sch,
        :netw_gq => netw_gq,
        :loc_matrix_res => loc_matrix_res,
        :loc_matrix_work => loc_matrix_work,
        :loc_lookup_res => loc_lookup_res,
        :loc_lookup_work => loc_lookup_work,
        :geo_lookup => geo_lookup,
        :k12_workers => k12_workers,
        :age_65o => age_65o,
        :mean_hh_connections => mean_hh_connections,
        :mean_wp_connections => mean_wp_connections,
        :distr_fn_loc_res => :exp, :distr_params_loc_res => (2,), 
        :distr_fn_loc_work => :exp, :distr_params_loc_work => (1,));

## precalc_vacc()
p_vacc_date_age = dser_path("precalc/p_vacc_date_age.jlse")
p_reduced_date_age = dser_path("precalc/p_reduced_date_age.jlse")
vac_wk_loc_age = dser_path("precalc/vac_wk_loc_age.jlse")
vreduced_wk_loc_age = dser_path("precalc/vreduced_wk_loc_age.jlse")
p_idxs_outside = dser_path("precalc/p_idxs_outside.jlse");
age_17u_county = dser_path("precalc/age_17u_county.jlse");
age_18_49_county = dser_path("precalc/age_18_49_county.jlse");
age_50_64_county = dser_path("precalc/age_50_64_county.jlse");
age_65o_county = dser_path("precalc/age_65o_county.jlse");
## assume outside commuters are 18-49
vac_outside_wk = round.(Int, p_vacc_date_age[:,2] * length(p_idxs_outside))
reduced_outside_wk = round.(Int, p_reduced_date_age[:,2] * length(p_idxs_outside))


nreps = 360
## store output
reps = Dict(i=>Dict() for i in 1:nreps)
## each run writes to a separate simunit object
## the objects share data (everything in "refs"), but never modify the shared data
## so multi-threading is fine
rlock = ReentrantLock();
Threads.@threads for i in 1:nreps

    params = Dict(
                :VE_inf => param_sets[i].VE,
                :pHosp => [param_sets[i].pH64u, param_sets[i].pH65o],
                :pResist => param_sets[i].S0, 
                :init_inf => param_sets[i].I0, 
                :pInf0 => 0.10, 
                :pInfAmp => 0.06, 
                :tPeak => 100)

    ## re-shuffle vaccinated each rep
    assign_17u = [shuffle(v) for v in age_17u_county]
    assign_18_49 = [shuffle(v) for v in age_18_49_county]
    assign_50_64 = [shuffle(v) for v in age_50_64_county]
    assign_65o = [shuffle(v) for v in age_65o_county]
    assign_outside = shuffle(p_idxs_outside);
    groups = [assign_17u, assign_18_49, assign_50_64, assign_65o];

    ## parallel scenarios with same params
    ## scenario A, normal vacc
    scen = :A
    vacSet_by_wk = assign_vac(groups, vac_wk_loc_age, assign_outside, vac_outside_wk)
    inputs = modelInputs(;vaccinated=vacSet_by_wk, params..., refs...)
    u = SimUnit(inputs)
    s = run(u,366)
    @lock rlock reps[i][scen] = s
    ## scenario B, reduced vacc
    scen = :B
    vacSet_by_wk = assign_vac(groups, vreduced_wk_loc_age, assign_outside, reduced_outside_wk)
    inputs = modelInputs(;vaccinated=vacSet_by_wk, params..., refs...)
    u = SimUnit(inputs)
    GC.gc()
    s = run(u,366)
    @lock rlock reps[i][scen] = s
    ## scenario C, no vacc
    scen = :C
    wks = keys(vacSet_by_wk)
    vacSet_by_wk = Dict(k=>Set(UInt32[]) for k in wks)
    inputs = modelInputs(;vaccinated=vacSet_by_wk, params..., refs...)
    u = SimUnit(inputs)
    GC.gc()
    s = run(u,366)
    @lock rlock reps[i][scen] = s

end

ser_path("precalc/reps.jlse", reps)





function hubformat(sc_idx)
    M = Float64.(round.(Int, matrices[sc_idx] * popsize / 100000))
    df = (DataFrame(M[1:hend,:],string.(1:nreps)))
    df[!,:horizon] = Vector{Union{Int,Missing}}(collect(1:hend))
    df = stack(df,1:nreps,"horizon")
    df[!,:scenario_id] .= scenario_ids[sc_idx]
    df[!,:target] .= "inc hosp"
    for i in 1:nreps
        push!(df, Dict(:horizon=>missing,:variable=>string(i),:value=>S0_rep[i],:scenario_id=>scenario_ids[sc_idx],:target=>"S0"))
    end
    df[!,:run_grouping] = parse.(Int, df[!,:variable])
    df[!,:stochastic_run] = length(scenario_ids)*(df[!,:run_grouping] .- 1) .+ sc_idx
    return df
end
