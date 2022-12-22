using Memento
const _LOGGER = Memento.getlogger(@__MODULE__)
var(pm::AbstractPowerModel, nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, nw)
var(pm::AbstractPowerModel, nw::Int, key::Symbol) = _IM.var(pm, pm_it_sym, nw, key)
var(pm::AbstractPowerModel, nw::Int, key::Symbol, idx) = _IM.var(pm, pm_it_sym, nw, key, idx)
var(pm::AbstractPowerModel, key::Symbol; nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, key; nw = nw)
var(pm::AbstractPowerModel, key::Symbol, idx; nw::Int=nw_id_default) = _IM.var(pm, pm_it_sym, key, idx; nw = nw)

#=
INFO:

- build_pf_all_flex_lazy takes into account the possibility to have downards and upwards flexibility. To do this, there is a linearization of the objective function
which would need the absolute value otherwise. The linearization is represented by the dummy variables and constraints. 

- this methods follows the same algorithm as in all the others _lazy mathods where congestion or voltage constraints are added iteratively if violated 

- here everything is included, so presence of DGs is tolerated (curtailment not yet added), upwards/downwards flexibility 

REMARKS:

- flexibility is introduced in the variable_load constraints.
- Flexibility and congestion coefficients are passed in the main code by using the functions add_load_flexibility and add_congestion_capacity

!!!! Here voltage is unbounded !!!!

=#

function build_pf_all_flex_lazy(pm::AbstractPowerModel)
    
    variable_bus_voltage(pm, bounded = false)
    variable_gen_power(pm, bounded = false)
    variable_branch_power(pm, bounded = false)
    variable_load(pm, bounded=true)
    variable_dummy(pm)

    objective_min_load_variations(pm)

    constraint_model_voltage(pm)  #do nothing

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3

        constraint_theta_ref(pm, i)
        constraint_voltage_magnitude_setpoint(pm, i)

    end

    for (i,bus) in ref(pm, :bus)
        constraint_power_balance_me(pm, i)

        # PQ Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            for j in ref(pm, :bus_gens, i)
                constraint_gen_setpoint_active(pm, j)
                constraint_gen_setpoint_reactive(pm, j)
            end
        end
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_to(pm,i)
        constraint_ohms_yt_from(pm,i)
    end

    for i in ids(pm, :load)
        constraint_load_factor(pm,i)
        constraint_dummy(pm,i) 
    end

end

function constraint_power_balance_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    #desa = ref(pm, nw, :load, ref(pm, nw, :bus_loads, i)[1], "pd")
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)
    bus_storage = ref(pm, nw, :bus_storage, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance_me_2(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
end

function constraint_power_balance_me_2(pm::AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
    vm   = var(pm, n, :vm, i)
    p    = get(var(pm, n),    :p, Dict()); PowerModels._check_var_keys(p, bus_arcs, "active power", "branch") #real power flowing from all branches (i to j)
    q    = get(var(pm, n),    :q, Dict()); PowerModels._check_var_keys(q, bus_arcs, "reactive power", "branch")
    
    load_p = get(var(pm, n), :load_p, Dict())
    load_q = get(var(pm, n), :load_q, Dict())
    
    pg   = get(var(pm, n),   :pg, Dict()); PowerModels._check_var_keys(pg, bus_gens, "active power", "generator")  #real power injected by generator
    qg   = get(var(pm, n),   :qg, Dict()); PowerModels._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, n),   :ps, Dict()); PowerModels._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n),   :qs, Dict()); PowerModels._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n),  :psw, Dict()); PowerModels._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n),  :qsw, Dict()); PowerModels._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(var(pm, n), :p_dc, Dict()); PowerModels._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(var(pm, n), :q_dc, Dict()); PowerModels._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    # the check "typeof(p[arc]) <: JuMP.NonlinearExpression" is required for the
    # case when p/q are nonlinear expressions instead of decision variables
    # once NLExpressions are first order in JuMP it should be possible to
    # remove this.
    nl_form = length(bus_arcs) > 0 && (typeof(p[iterate(bus_arcs)[1]]) <: JuMP.NonlinearExpression)

    if !nl_form
        cstr_p = JuMP.@constraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(load_p[l] for l in bus_loads)
            - sum(gs for (i,gs) in bus_gs)*vm^2
        )
    else
        cstr_p = JuMP.@NLconstraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(load_p[l] for l in bus_loads)
            - sum(gs for (i,gs) in bus_gs)*vm^2
        )
    end

    if !nl_form
        cstr_q = JuMP.@constraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(load_q[l] for l in bus_loads)
            + sum(bs for (i,bs) in bus_bs)*vm^2
        )
    else
        cstr_q = JuMP.@NLconstraint(pm.model,
            sum(q[a] for a in bus_arcs)
            + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
            ==
            sum(qg[g] for g in bus_gens)
            - sum(qs[s] for s in bus_storage)
            - sum(load_q[l] for l in bus_loads)
            + sum(bs for (i,bs) in bus_bs)*vm^2
        )
    end

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

function variable_load(pm::AbstractACPModel; kwargs...)
    variable_load_real(pm; kwargs...)
    variable_load_imm(pm; kwargs...)
end

function variable_load_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

    load_p= var(pm,nw)[:load_p] = JuMP.@variable(pm.model, y[i in ids(pm, nw, :load)],
    base_name="$(nw)_load_p",  start = comp_start_value(ref(pm, nw, :load, i), "pd"))

    if bounded
        for (i, load) in ref(pm, nw, :load)
            JuMP.set_lower_bound(load_p[i], (1-load["flex_%"]) * load["pd"]) # felxibility 
            JuMP.set_upper_bound(load_p[i], (1+load["flex_%"]) * load["pd"])
        end
    end

    report && sol_component_value(pm, nw, :load, :load_p, ids(pm, nw, :load), load_p)

end

function variable_load_imm(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    
    load_q = var(pm,nw)[:load_q] = JuMP.@variable(pm.model, x[i in ids(pm, nw, :load)],
    base_name="$(nw)_load_q",  start = comp_start_value(ref(pm, nw, :load, i), "qd"))

    if bounded
        for (i, load) in ref(pm, nw, :load)
            JuMP.set_lower_bound(load_q[i], (1-load["flex_%"]) * load["qd"]) # felxibility 
            JuMP.set_upper_bound(load_q[i], (1+load["flex_%"]) * load["qd"])
        end
    end

    report && sol_component_value(pm, nw, :load, :load_q, ids(pm, nw, :load), load_q)

end

function objective_min_load_variations(pm::AbstractPowerModel; nw = nw_id_default, kwargs...)
    
    return JuMP.@objective(pm.model, Min, sum(var(pm,nw, :x_p, i) + var(pm,nw, :y_p, i) for i in ids(pm, :load)) + sum(var(pm,nw, :x_q, i) + var(pm,nw, :y_q, i) for i in ids(pm, :load)))
end

function constraint_thermal_limit_from_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        p_fr = var(pm, nw, :p, f_idx)
        q_fr = var(pm, nw, :q, f_idx)
        
        JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= (branch["cong_cap"]*branch["rate_a"])^2)  #80% loading 
    end
end

function constraint_thermal_limit_to_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if haskey(branch, "rate_a")
        p_to = var(pm, nw, :p, t_idx)
        q_to = var(pm, nw, :q, t_idx)

        JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= (branch["cong_cap"] * branch["rate_a"])^2)  #80% loading 
    end
end

function constraint_voltage_magnitude_upper_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    v_max = bus["vmax"]

    v = var(pm, nw, :vm, i)
    JuMP.@constraint(pm.model, v <= v_max)

end

function constraint_voltage_magnitude_lower_me(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    v_min = bus["vmin"]

    v = var(pm, nw, :vm, i)
    JuMP.@constraint(pm.model, v_min <= v)

end

function variable_dummy(pm::AbstractPowerModel, nw::Int=nw_id_default, report::Bool=true)

    var(pm,nw)[:y_p] = JuMP.@variable(pm.model, 0 <= y_p[i in ids(pm, nw, :load)])
    var(pm,nw)[:x_p] = JuMP.@variable(pm.model, 0 <= x_p[i in ids(pm, nw, :load)])
    var(pm,nw)[:y_q] = JuMP.@variable(pm.model, 0 <= y_q[i in ids(pm, nw, :load)])
    var(pm,nw)[:x_q] = JuMP.@variable(pm.model, 0 <= x_q[i in ids(pm, nw, :load)])

    report && sol_component_value(pm, nw, :load, :x_p, ids(pm, nw, :load), x_p)
    report && sol_component_value(pm, nw, :load, :x_q, ids(pm, nw, :load), x_q)
    report && sol_component_value(pm, nw, :load, :y_p, ids(pm, nw, :load), y_p)
    report && sol_component_value(pm, nw, :load, :y_q, ids(pm, nw, :load), y_q)

end

function constraint_load_factor(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    load = ref(pm,nw,:load,i)
    load_p = get(var(pm, nw), :load_p, Dict())
    load_q = get(var(pm, nw), :load_q, Dict())

    cosφ = load["pd"]/sqrt(load["pd"]^2+load["qd"]^2)
    φ = acos(cosφ)
    tanφ = tan(φ)

    pd = load_p[i]
    qd = load_q[i]

    @constraint(pm.model, pd*tanφ == qd)
end


function constraint_dummy(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
    load = ref(pm,nw,:load,i)
    load_p = get(var(pm, nw), :load_p, Dict())
    load_q = get(var(pm, nw), :load_q, Dict())
    y_p = get(var(pm,nw), :y_p, Dict())
    x_p = get(var(pm,nw), :x_p, Dict())
    y_q = get(var(pm,nw), :y_q, Dict())
    x_q = get(var(pm,nw), :x_q, Dict())

    p = load["pd"]
    p_flex = load_p[i]
    JuMP.@constraint(pm.model, p - p_flex <= y_p[i])
    JuMP.@constraint(pm.model, p_flex - p <= x_p[i])

    q = load["qd"]
    q_flex = load_q[i]
    JuMP.@constraint(pm.model, q - q_flex <= y_q[i])
    JuMP.@constraint(pm.model, q_flex - q <= x_q[i])
end