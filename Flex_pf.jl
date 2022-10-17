using PowerModels
using PowerModelsAnalytics
using PowerPlots
using ColorSchemes
using JuMP, Ipopt
using Setfield
using Plots
import InfrastructureModels; const _IM = InfrastructureModels

#=

    Compute PF by using either: my_build_pf or my_build_pf_DGs. The former takes generators as PV type. The latter as PQ. 

=#


include("Network_handling.jl")
include("Network_paths.jl")
include("My_ref//My_ref.jl")
include("My_functions.jl")

# Input file
file_name = "Official_rural.m"
file_path = "C://Workdir//Develop//Official_rural.m"
net_data = parse_file(file_path)

# Add flexibility % that each load can offer
flex = add_load_flexibility(net_data, 0.2)
net_data["flex"] = 20
update_data!(net_data, flex)

# Add congestion capacity for each line
cong_cap = add_congestion_capacity(net_data, 1.0)
update_data!(net_data, cong_cap)

# Solve model
pm = instantiate_model(net_data, ACPPowerModel, my_build_pf)
result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))

# Compute flexibility offered by each load 
flex_loads, p_load, q_load = calc_flexibility_offered(net_data, result)

# Compute branch flows and losses
update_data!(net_data, result["solution"])
flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)

losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)

# Voltage profile with flexibility 

feeder_ID, path_volt = calc_voltage_profile(net_data, result, file_name)

# Print statements

print("\n Real power generated: ", result["solution"]["gen"]["1"]["pg"], "MW \n")
print("\n Immaginary power generated: ", result["solution"]["gen"]["1"]["qg"], "MVar \n")

print("\n Real power requested: ", p_load, "MW \n")
print("\n Immaginary power requested: ", q_load, "MVar \n")

print("\n Real power losses: ", sum(p_loss), "MW \n")
print("\n Immaginary power losses: ", sum(q_loss), "MVar \n")

print("\n Flexibility offered: ", sum([flex["diff_real"] for (i,flex) in flex_loads]), "MW \n")

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

# Look at finction description in My_functions to see which argument you can pass
plot_grid(net_data, "flex")


