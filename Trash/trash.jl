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

include("My_ref/My_ref_no_flex_with_DG.jl")
include("C:/Workdir/Develop/PF_simulations/My_ref//My_ref_lazy_flex.jl")  
include("C:/Workdir/Develop/PF_simulations/My_functions.jl")
include("C:/Workdir/Develop/Hosting_Capacity/HC_functions.jl")


#Parameters

congestion_limit = 100  	             # congestion limit %
threshold = 100                        # value used to identify branches with current rating higher than threshold


size = 11                              # size of generator/s
seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 5              # number of random DGs per feeder
power_target_per_feeder = 10            # total capacity installed in each feeder

# Input file
file_name = "Official_urban.m"
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
pm = instantiate_model(net_data, ACPPowerModel, no_flex_with_DG)
result = optimize_model!(pm, optimizer=Ipopt.Optimizer)
update_data!(net_data, result["solution"])

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

plot_grid(net_data, "basic", "basic","loading")