using PowerModels
using PowerModelsAnalytics
using PowerPlots
using ColorSchemes
using JuMP, Ipopt
using Setfield
using Plots
import InfrastructureModels; const _IM = InfrastructureModels

include("My_ref//test_ref.jl")
include("My_ref//My_ref_lazy_DGs.jl")
include("My_ref//My_ref_DGs.jl")
include("My_functions.jl")

# Input file
file_name= "Official_rural.m"
file_path = "C://Users//u0152683//Desktop//Networks//PF simulation//Official_rural.m"
net_data = parse_file(file_path)

# Add flexibility % that each load can offer
flex = add_load_flexibility(net_data, 0.2)
net_data["flex"] = 20
update_data!(net_data, flex)

# Add congestion capacity for each line
cong_cap = add_congestion_capacity(net_data, 0.8)
update_data!(net_data, cong_cap)


# Solve
result = solve_pf_branch_power_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_all_flex_lazy)
#pm = instantiate_model(net_data, ACPPowerModel, my_build_pf_DGs)
#result = optimize_model!(pm, optimizer=Ipopt.Optimizer)
#ciao = solve_pf_branch_power_cuts_mine(net_data, ACPPowerModel, Ipopt.Optimizer, build_pf_lazy_DGs)

# Compute flexibility offered by each load 
flex_loads, p_load, q_load = flexibility_offered_new(net_data, result)

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

