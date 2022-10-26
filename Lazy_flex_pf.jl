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

#include("My_ref//My_ref_lazy_DGs.jl")
include("My_ref//My_ref_lazy_flex.jl")  
include("My_functions.jl")
include("C:/Workdir/Develop/Hosting_Capacity/HC_function_curt.jl")

#=

    Compute PF by using either: 
    -   build_pf_lazy: takes generators as PV. Only downwards flexibility. NO DG curtailment
    -   build_pf_lazy_DGs: takes generators as PQ. Only downwards flexibility. NO DG curtailment
    -   build_pf_all_flex_lazy: takes generators as PQ. Upwards and downwards flexibility included. NO DG curtailment 

=#


#Parameters

flex = 20                              # flexibility offered %
congestion_limit = 100  	             # congestion limit %
threshold = 100                        # value used to identify branches with current rating higher than threshold
curtailment = 0                       # DG curtailment %

size = 11                              # size of generator/s
seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 5              # number of random DGs per feeder
power_target_per_feeder = 10            # total capacity installed in each feeder

# Input file
file_name = "Official_rural.m"
file_path = "C://Workdir//Develop//Official_rural.m"
net_data = parse_file(file_path)

# Add flexibility % that each load can offer
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

size_std = power_target_per_feeder/gen_number_per_feeder  #size of each DG
add_generators(net_data, random_generators, size_std, curtailment/100)

#= Add a random generator
Random.seed!(seed)
passive_nodes = collect(keys(net_data["bus"]))
deleteat!(passive_nodes, findall(x->x=="1", passive_nodes))
deleteat!(passive_nodes, findall(x->x=="$mv_busbar", passive_nodes))
active_node = sample(passive_nodes,1)[1]    # if more than one samples are done, then remove the [1] 

add_single_generator(net_data, size, active_node, curtailment/100)
=#
gen_ID = get_gen_info(net_data, feeder_ID_1)

# Solve PF 
#result = solve_pf_branch_power_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_all_flex_lazy)
result = solve_pf_branch_voltage_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_all_flex_lazy)
update_data!(net_data, result["solution"])

# Compute flexibility offered by each load 
flex_loads, p_load, q_load = calc_flexibility_offered_new(net_data, result)
#flex_loads, p_load, q_load = calc_flexibility_offered(net_data, result)

# Compute branch flows and losses

flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)

losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)

# Voltage profile with flexibility 

feeder_ID, path_volt, mv_busbar = calc_voltage_profile(net_data, result, file_name)
gen_ID = get_gen_info(net_data, feeder_ID)

# Update net_data with branch loadings 
calc_branch_loading(net_data, feeder_ID, gen_ID, threshold)

#Printing statements
printing_statements(result, p_load, q_load, p_loss, q_loss, flex_loads, gen_ID, flows, threshold, mv_busbar)

# Create df for storing data
df = create_df(feeder_ID, result, threshold)
print(df)

#Look at function description in My_functions to see which argument you can pass
#bus, gen, branch
plot_grid(net_data, "p_flex","basic", "loading"; display_flow = true)
