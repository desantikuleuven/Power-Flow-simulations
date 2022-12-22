using PowerModels
using PowerModelsAnalytics
using PowerPlots
using ColorSchemes
using JuMP, Ipopt
using Setfield
using Plots
using VegaLite
import InfrastructureModels; const _IM = InfrastructureModels
using DataFrames
using StatsBase
using Random

#=
    Script to perform PF without any kind of flexibility.
    If solved by using no_flex_with_DG, the DGs are considered as PQ buses with VOLTAGE CONSTRAINTS. 
    Otherwise use simply solve_ac_pf
=#

include("My_ref/My_ref_no_flex_with_DG.jl")
include("My_functions.jl")
include("C:/Workdir/Develop/Hosting_Capacity/HC_functions.jl")


#Parameters

congestion_limit = 100  	             # congestion limit %
threshold = 100                        # value used to identify branches with current rating higher than threshold


size = 11                              # size of generator/s
seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 5              # number of random DGs per feeder
power_target_per_feeder = 10            # total capacity installed in each feeder

# Input file
file_name = "Official_rural.m"
file_path = "C://Workdir//Develop//"*file_name
net_data = parse_file(file_path)

# Add congestion capacity for each line

cong_cap = add_congestion_capacity(net_data, congestion_limit/100)
update_data!(net_data, cong_cap)

# Get feeder info 
feeder_ID_1, mv_busbar = get_feeder_data(net_data, file_name)

#Random choice of buses
random_generators = get_random_generators(feeder_ID_1, gen_number_per_feeder, seed)  

# Add new generators
#add_generators(net_data, random_generators, size)

# Get random generator
active_node = get_random_DG(seed,net_data)

# Add a random generator
add_single_generator(net_data, size, active_node)

# Create dict for DGs (CALL ONLY BEFORE PF)
gen_ID = get_gen_info(net_data, feeder_ID_1)

# Run PF

model = "no_flex_with_DG"
# result = solve_ac_pf(net_data, Ipopt.Optimizer)
pm = instantiate_model(net_data, ACPPowerModel, no_flex_with_DG)
result = optimize_model!(pm, optimizer=Ipopt.Optimizer)
update_data!(net_data, result["solution"])

@assert result["termination_status"] == LOCALLY_SOLVED

# Compute branch flows and losses
flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)

losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)

# Voltage profile with flexibility 
save_path =  create_save_path(file_name, gen_ID, model)
feeder_ID, path_volt, mv_busbar = calc_voltage_profile(net_data, result, file_name, save_path, false,true,true)

# Update net_data with branch loadings 
calc_branch_loading(net_data, feeder_ID, gen_ID, threshold)


p_load = 0
q_load = 0

for (i,load) in net_data["load"]
    p_load += load["pd"]
    q_load += load["qd"]
end


p = []
q = []
for (ID,data) in flows["branch"]
    push!(p,data["pt"])
    push!(q,data["qt"])
end

print("\n Real power generated: ", result["solution"]["gen"]["1"]["pg"], "MW \n")
print("\n Immaginary power generated: ", result["solution"]["gen"]["1"]["qg"], "MVar \n")

print("\n Real power requested: ", p_load, "MW \n")
print("\n Immaginary power requested: ", q_load, "MVar \n")

print("\n Real power losses: ", sum(p_loss), "MW \n")
print("\n Immaginary power losses: ", sum(q_loss), "MVar \n")

ratings = []
[push!(ratings, data["rate_a"]) for (ID,data) in net_data["branch"]]

s = abs.(complex.(p,q))  # Absolute value of Apparent power flowing in each branch
branch_loading = round.(s./ratings, digits = 3)*100
index = 80
num = length(filter(x->x>index, branch_loading))

println("\nBranches with loading higher than ", index, "%: ", num)


plot_grid(net_data, "basic", "basic","loading")