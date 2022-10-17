using PowerModels
using PowerModelsAnalytics
using PowerPlots
using ColorSchemes
using JuMP, Ipopt
using Setfield
using Plots
import InfrastructureModels; const _IM = InfrastructureModels
using DataFrames

#include("My_ref//My_ref_lazy_DGs.jl")
include("My_ref//My_ref_lazy_ALL.jl")  
include("My_functions.jl")

#=

    Compute PF by using either: 
    -   build_pf_lazy: takes generators as PV. Only downwards flexibility. NO DG curtailment
    -   build_pf_lazy_DGs: takes generators as PQ. Only downwards flexibility. NO DG curtailment
    -   build_pf_all_flex_lazy: takes generators as PQ. Upwards and downwards flexibility included. NO DG curtailment 

=#

# Input file
file_name = "Official_rural_test.m"
file_path = "C://Users//u0152683//Desktop//Networks//PF simulation//Official_rural_test.m"
net_data = parse_file(file_path)

# Add flexibility % that each load can offer
net_data["flex"] = 20
flex = add_load_flexibility(net_data, net_data["flex"]/100)
update_data!(net_data, flex)

# Add congestion capacity for each line
congestion_limit = 300
threshold = 100
cong_cap = add_congestion_capacity(net_data, congestion_limit/100)
update_data!(net_data, cong_cap)

# Solve PF 
result = solve_pf_branch_power_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_all_flex_lazy)
#result = solve_pf_branch_voltage_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_all_flex_lazy)
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

feeder_ID, path_volt, mv_busbar = calc_voltage_profile(net_data, result, file_name, false,false)
gen_ID = get_gen_info(net_data, feeder_ID)

# Update net_data with branch loadings 
calc_branch_loading(net_data, feeder_ID, gen_ID, threshold)

#Printing statements
printing_statements(result, p_load, q_load, p_loss, q_loss, flex_loads, up_p, up_q, down_p, down_q, gen_ID, flows, threshold, mv_busbar)

# Create df for storing data
df = create_df(feeder_ID, result, threshold)
print(df)

#Look at function description in My_functions to see which argument you can pass
#plot_grid_new(net_data, "basic")

