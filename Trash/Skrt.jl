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

include("My_ref/My_ref_lazy_DG_curtail.jl")
include("My_functions.jl")
include("C:/Workdir/Develop/Hosting_Capacity/HC_function_curt.jl")

#Parameters

flex = 20                              # flexibility offered %
congestion_limit = 100  	             # congestion limit %
threshold = 100                        # value used to identify branches with current rating higher than threshold
curtailment = 20                       # DG curtailment %

size = 2                              # size of generator/s
seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 5             # number of random DGs per feeder
power_target_per_feeder = 5            # total capacity installed in each feeder


# Input file
file_name = "Official_urban.m"
file_path = "C://Workdir//Develop//"*file_name
net_data = parse_file(file_path)

# Add flexibility % that each load can offer
net_data["curtailment"] = curtailment
net_data["flex"] = flex
update_data!(net_data, add_load_flexibility(net_data, net_data["flex"]/100))

# Add congestion capacity for each line

cong_cap = add_congestion_capacity(net_data, congestion_limit/100)
update_data!(net_data, cong_cap)

# Get feeder info 
feeder_ID_1, mv_busbar = get_feeder_data(net_data, file_name)

#Random choice of buses
random_generators = get_random_generators(feeder_ID_1, gen_number_per_feeder, seed)  

# Add new generators
#add_generators(net_data, random_generators, size, curtailment/100)

# Get random generator
active_node = get_random_DG(seed,net_data)

# Add a random generator
add_single_generator(net_data, size, active_node; curt =  curtailment/100)

# Create dict for DGs (CALL ONLY BEFORE PF)
gen_ID = get_gen_info(net_data, feeder_ID_1)

# Solve PF 
model = "build_pf_dg_curtail" 
#result = solve_pf_branch_power_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_dg_curtail)
result = solve_pf_branch_voltage_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_dg_curtail)
update_data!(net_data, result["solution"])

# Compute flexibility offered by each load 
flex_loads, p_load, q_load = calc_flexibility_offered_final(net_data, result)

# Compute curtailment for each generator, update net_data and gen_ID
feeder_curtailment = feeder_dg_curtailment(net_data, result, gen_ID)

# Compute branch flows and losses
flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)

losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)

# Voltage profile with flexibility 
save_path =  create_save_path(file_name, gen_ID, model)
feeder_ID, path_volt, mv_busbar = calc_voltage_profile(net_data, result, file_name, save_path, true, true, true)

# Update net_data with branch loadings 
calc_branch_loading(net_data, feeder_ID, gen_ID, threshold)

#Printing statements
printing_statements(result, p_load, q_load, p_loss, q_loss, flex_loads, gen_ID, flows, threshold, mv_busbar)

#Look at function description in My_functions to see which argument you can pass
#bus, gen, branch
plot_grid(net_data, "basic","basic", "loading"; display_flow = false)




