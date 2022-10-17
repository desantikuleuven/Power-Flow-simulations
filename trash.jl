using PowerModels
using PowerModelsAnalytics
using PowerPlots
using ColorSchemes
using JuMP, Ipopt
using Setfield
using Plots
import InfrastructureModels; const _IM = InfrastructureModels

#=
6831 4 0.0 0.0 0.0 1.0 1.0 1 7.0 0.0 0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	
6807 4 0.0 0.0 0.0 1.0 1.0 1 7.0 0.0 0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
6607 4 0.0 0.0 0.0 1.0 1.0 1 7.0 0.0 0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
6867 4 0.0 0.0 0.0 1.0 1.0 1 7.0 0.0 0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
=#

include("My_ref//test_ref.jl")
include("My_functions.jl")

# Input file
file_name = "Official_rural_test.m"
file_path = "C://Users//u0152683//Desktop//Networks//PF simulation//Official_rural_test.m"
net_data = parse_file(file_path)

# Add flexibility % that each load can offer
net_data["flex"] = 0
flex = add_load_flexibility(net_data, net_data["flex"]/100)
update_data!(net_data, flex)

# Add congestion capacity for each line
cong_cap = add_congestion_capacity(net_data, 5.0)
update_data!(net_data, cong_cap)


# Solve PF 
result = solve_pf_branch_power_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_all_flex_lazy)
#result = solve_pf_branch_voltage_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_all_flex_lazy)

# Compute flexibility offered by each load 
flex_loads, p_load, q_load = calc_flexibility_offered_new(net_data, result)

# Compute branch flows and losses
update_data!(net_data, result["solution"])
flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)


losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)

# Voltage profile with flexibility 
feeder_ID, path_volt = calc_voltage_profile(net_data, result, file_name)

gen_ID = get_gen_info(net_data, feeder_ID)

# Update net_data with branch loadings 
calc_branch_loading(net_data, feeder_ID,gen_ID)


# Print statements
print("\n Termination status: ", result["termination_status"])
print("\n Objective function: ", result["objective"])

print("\n Real power generated: ", result["solution"]["gen"]["1"]["pg"], "MW \n")
print("\n Immaginary power generated: ", result["solution"]["gen"]["1"]["qg"], "MVar \n")

print("\n Real power requested: ", p_load, "MW \n")
print("\n Immaginary power requested: ", q_load, "MVar \n")

print("\n Real power losses: ", sum(p_loss), "MW \n")
print("\n Immaginary power losses: ", sum(q_loss), "MVar \n")

print("\n Active power Flexibility offered: ", sum([abs(flex["diff_real"]) for (i,flex) in flex_loads]), "MW \n")
print("\n Rective power Flexibility offered: ", sum([abs(flex["diff_imm"]) for (i,flex) in flex_loads]), "MVar \n")

up_p = [flex["diff_real"] for (i,flex) in flex_loads if flex["diff_real"]>0]
up_q = [flex["diff_imm"] for (i,flex) in flex_loads if flex["diff_imm"]>0]

down_p = [flex["diff_real"] for (i,flex) in flex_loads if flex["diff_real"]<0]
down_q = [flex["diff_imm"] for (i,flex) in flex_loads if flex["diff_imm"]<0]

if length(up_p)>0 
    print("\n P Upwards Flexibility offered: ", sum(up_p), "MW \n")
end
if length(down_p) > 0 
    print("\n P Downards Flexibility offered: ", sum(down_p), "MW \n")
end
if length(up_q)>0 
    print("\n Q Upwards Flexibility offered: ", sum(up_q), "MVar \n")
end
if length(down_q) > 0 
    print("\n Q Downards Flexibility offered: ", sum(down_q), "MVar \n")
end

[println("\n Vmin: ", round(feeder_ID[gen["ref"]]["vmin"], digits = 5)," Vmax: ",round(feeder_ID[gen["ref"]]["vmax"],digits = 5)," for ",gen["feeder"]," (where generator $i is connected)") for (i,gen) in gen_ID]
[println("\n Feeder loading: ", round(gen["max_branch_load"],digits = 3),"% in branch ",gen["Critical_branch"]) for (i,gen) in gen_ID]

#Branch loading info
p = []
q = []
for (ID,data) in flows["branch"]
    push!(p,data["pt"])
    push!(q,data["qt"])
end

ratings = []
[push!(ratings, data["rate_a"]) for (ID,data) in net_data["branch"]]

s = abs.(complex.(p,q))  # Absolute value of Apparent power flowing in each branch
branch_loading = round.(s./ratings, digits = 3)*100
index = 70
num = length(filter(x->x>index, branch_loading))

println("\nBranches with loading higher than ", index, "%: ", num)

vmin = Vector()
vmax = Vector()
[push!(vmin,feeder["vmin"]) for (id, feeder) in feeder_ID ]
println("\n Lowest voltage magnitude in the grid: ", minimum(vmin))
[push!(vmin,feeder["vmax"]) for (id, feeder) in feeder_ID ]
println("\n Highest voltage magnitude in the grid: ", maximum(vmin))

plot_grid_new(net_data, "basic", true)
#=
function sensitivity(net_data::Dict, P::Float64,i)

    println("\n\n")
    println("********** ITERATION $i: P= $P MW**************")
    println("\n\n")

    net_data["gen"]["2"]["pg"] = P

    # Solve PF 
    result = solve_pf_branch_power_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_lazy_DGs)
    #result = solve_pf_branch_voltage_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_lazy_DGs)

    # Compute flexibility offered by each load 
    flex_loads, p_load, q_load = flexibility_offered(net_data, result)

    # Compute branch flows and losses
    update_data!(net_data, result["solution"])
    flows = calc_branch_flow_ac(net_data)
    update_data!(net_data, flows)


    losses, p_loss, q_loss = calc_power_losses(net_data)
    update_data!(net_data, losses)

    # Voltage profile with flexibility 
    feeder_ID, path_volt = calc_voltage_profile(net_data, result, file_name)

    gen_ID = get_gen_info(net_data, feeder_ID)

    # Update net_data with branch loadings 
    calc_branch_loading(net_data, feeder_ID,gen_ID)


    # Print statements

    println("\n MV busbar voltage: ", net_data["bus"]["6793"]["vm"])

    print("\n Flexibility offered: ", sum([flex["diff_real"] for (i,flex) in flex_loads]), "MW \n")


    [println("\n Vmin: ", round(feeder_ID[gen["ref"]]["vmin"], digits = 5)," Vmax: ",round(feeder_ID[gen["ref"]]["vmax"],digits = 5)," for ",gen["feeder"]," (where generator $i is connected)") for (i,gen) in gen_ID]
    [println("\n Feeder loading: ", round(gen["max_branch_load"],digits = 3),"% in branch ",gen["Critical_branch"]) for (i,gen) in gen_ID]

    plot_grid(net_data, "flex")
end

for (i, P) in enumerate([10.0, 14.0])
    sensitivity(net_data, P, i)
end
=#
