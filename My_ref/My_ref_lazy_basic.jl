using Memento
const _LOGGER = Memento.getlogger(@__MODULE__)

function build_pf_lazy(pm::AbstractPowerModel)
    
    variable_bus_voltage(pm, bounded = false)
    variable_gen_power(pm, bounded = false)
    variable_branch_power(pm, bounded = false)
    variable_load(pm, bounded=true)

    objective_min_load_variations(pm)

    constraint_model_voltage(pm)  #do nothing

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3

        constraint_theta_ref(pm, i)
        constraint_voltage_magnitude_setpoint(pm, i)

    end

    for (i,bus) in ref(pm, :bus)
        constraint_power_balance_me(pm, i)

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_voltage_magnitude_setpoint(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_gen_setpoint_active(pm, j)
            end
        end
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_to(pm,i)
        constraint_ohms_yt_from(pm,i)
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

    load_p = var(pm,nw)[:load_p] = JuMP.@variable(pm.model, [i in ids(pm, nw, :load)],
    base_name="$(nw)_load_p",  start = comp_start_value(ref(pm, nw, :load, i), "pd"))

    if bounded
        for (i, load) in ref(pm, nw, :load)
            JuMP.set_lower_bound(load_p[i], (1-load["flex_%"]) * load["pd"]) # felxibility 
            JuMP.set_upper_bound(load_p[i], load["pd"])
        end
    end

    report && sol_component_value(pm, nw, :load, :load_p, ids(pm, nw, :load), load_p)

end

function variable_load_imm(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    
    load_q = var(pm,nw)[:load_q] = JuMP.@variable(pm.model, [i in ids(pm, nw, :load)],
    base_name="$(nw)_load_q",  start = comp_start_value(ref(pm, nw, :load, i), "qd"))

    if bounded
        for (i, load) in ref(pm, nw, :load)
            JuMP.set_lower_bound(load_q[i], (1-load["flex_%"]) * load["qd"]) # felxibility 
            JuMP.set_upper_bound(load_q[i], load["qd"])
        end
    end

    report && sol_component_value(pm, nw, :load, :load_q, ids(pm, nw, :load), load_q)

end

function objective_min_load_variations(pm::AbstractPowerModel; nw = nw_id_default, kwargs...)
    
    return JuMP.@objective(pm.model, Min, sum(load["pd"] - var(pm,nw, :load_p, i) for (i,load) in ref(pm, :load)) + sum(load["qd"] - var(pm,nw, :load_q, i) for (i,load) in ref(pm, :load)) )

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

#=

function vaiable_bus_voltage_me(pm::AbstractACPModel; kwargs...)
    variable_bus_voltage_angle_me(pm; kwargs...)
    variable_bus_voltage_magnitude_me(pm; kwargs...)
end

function variable_bus_voltage_angle_me(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    va = var(pm, nw)[:va] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_va",
        start = comp_start_value(ref(pm, nw, :bus, i), "va_start")
    )

    report && sol_component_value(pm, nw, :bus, :va, ids(pm, nw, :bus), va)
end

function variable_bus_voltage_magnitude_me(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vm = var(pm, nw)[:vm] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_vm",
        start = comp_start_value(ref(pm, nw, :bus, i), "vm_start", 1.0)
    )

    if bounded
        for (i, bus) in ref(pm, nw, :bus)
            JuMP.set_lower_bound(vm[i], bus["vmin"])
            JuMP.set_upper_bound(vm[i], bus["vmax"])
        end
    end

    report && sol_component_value(pm, nw, :bus, :vm, ids(pm, nw, :bus), vm)
end
=#