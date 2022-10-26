using PowerModels
using PowerModelsAnalytics
using PowerPlots
using PowerPlots.Experimental
using ColorSchemes
using JuMP, Ipopt
using Setfield
using Plots
using VegaLite
import InfrastructureModels; const _IM = InfrastructureModels
using DataFrames
using StatsBase
using Random

#include("Network_handling.jl")
#include("Network_paths.jl")
#include("My_ref.jl")

include("My_functions.jl")

file_name = "Official_rural.m"
file_path = "C://Workdir//Develop//Official_rural.m"
net_data = parse_file(file_path)

result = solve_ac_pf(net_data, Ipopt.Optimizer)
update_data!(net_data, result["solution"])

flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)

losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)

feeder_ID, path_volt = calc_voltage_profile(net_data, result, file_name, false)


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

#plot1 = powerplot(net_data; gen_size = 50, bus_size = 50, load_size = 10, branch_size = 3, width = 1200, height = 1200)

#bus, gen, branch 
plot_grid_now(net_data, "vm","curtailment","basic"; zoom =true, display_flow = false)

