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
include("My_ref//My_ref_DGs.jl")
include("My_functions.jl")
#Parameters

flex = 20                              # flexibility offered %
congestion_limit = 100  	             # congestion limit %
threshold = 100                        # value used to identify branches with current rating higher than threshold

size = 11                              # size of generator/s
seed = 99                              # seed for random DG buses choice
gen_number_per_feeder = 5              # number of random DGs per feeder
power_target_per_feeder = 10            # total capacity installed in each feeder

# Input file
# Input file
file_name = "Official_urban.m"
file_path = "C://Workdir//Develop//"*file_name
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
#add_generators(net_data, random_generators, size, curtailment/100)

# Get random generator
active_node = get_random_DG(seed,net_data)

# Add a random generator
add_single_generator(net_data, size, active_node)

# Create dict for DGs (CALL ONLY BEFORE PF)
gen_ID = get_gen_info(net_data, feeder_ID_1)

# Solve model
pm = instantiate_model(net_data, ACPPowerModel, my_build_pf_DGs)
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
save_path =  create_save_path(file_name, gen_ID, model)
feeder_ID, path_volt, mv_busbar = calc_voltage_profile(net_data, result, file_name, save_path)

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

#Look at function description in My_functions to see which argument you can pass
#bus, gen, branch
plot_grid(net_data, "p_flex","basic", "loading"; display_flow = true)


