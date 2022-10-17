using PowerModels
using DataFrames
using JuMP, Ipopt
import InfrastructureModels; const _IM = InfrastructureModels

include("Network_handling.jl")
include("Network_paths.jl")
#include("My_ref.jl")
include("My_functions.jl")

file_name = "Official_rural_test_dg.m"
file_path = "C://Users//u0152683//Desktop//Networks//PF simulation//Official_rural_test_dg.m"
net_data = parse_file(file_name)

function alternative_pf(pm::AbstractPowerModel)

    variable_bus_voltage(pm, bounded = false)
    variable_gen_power(pm, bounded = false)
    variable_branch_power(pm, bounded = false)

    constraint_model_voltage(pm)

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3

        constraint_theta_ref(pm, i)
        constraint_voltage_magnitude_setpoint(pm, i)

    end

    for (i,bus) in ref(pm, :bus)
        constraint_power_balance(pm, i)

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_voltage_magnitude_setpoint(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_gen_setpoint_active(pm, j)
            end
        end
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_to(pm,i)
        constraint_ohms_yt_from(pm,i)
    end
end

pm = instantiate_model(net_data, ACPPowerModel, alternative_pf)
result = optimize_model!(pm, optimizer=Ipopt.Optimizer)

update_data!(net_data, result["solution"])
flows = calc_branch_flow_ac(net_data)
update_data!(net_data, flows)

losses, p_loss, q_loss = calc_power_losses(net_data)
update_data!(net_data, losses)

p_load = 0
q_load = 0

for (i,load) in net_data["load"]
    p_load += load["pd"]
    q_load += load["qd"]
end

print("\n Real power generated: ", result["solution"]["gen"]["1"]["pg"], "MW \n")
print("\n Immaginary power generated: ", result["solution"]["gen"]["1"]["qg"], "MVar \n")

print("\n Real power requested: ", p_load, "MW \n")
print("\n Immaginary power requested: ", q_load, "MVar \n")

print("\n Real power losses: ", sum(p_loss), "MW \n")
print("\n Immaginary power losses: ", sum(q_loss), "MVar \n")

function plot_grid(net_data::Dict, bus_data_input::String, zoom = true)

    [net_data["branch"][i]["loading"] = round(abs(complex(branch["pt"],branch["qt"]))/branch["rate_a"], digits = 3)*100 for (i,branch) in net_data["branch"]]

    plot1 = powerplot(
        net_data,
        width = 750, 
        height = 750,
        # bus_data=:vm,
        # bus_data_type=:quantitative,
        gen_data=:pg,
        gen_data_type=:quantitative,
        branch_data = "loading",
        branch_data_type=:quantitative,
        bus_data = bus_data_input,
        bus_data_type = "quantitative",
        bus_color = ["#C0C0C0","#000000"],
        branch_color=["green","orange","red"],
        gen_color=["blue","purple"],
        #gen_color = "purple",
        branch_size = 2, 
        bus_size = 50,
        gen_size = 80,
    )

    
    plot1.layer[4]["transform"] = Dict{String, Any}[
    Dict("calculate"=>"datum.pg/datum.pmax*100", "as"=>"gen_Percent_Loading"),
    Dict("calculate"=>"datum.pg", "as"=>"GenPower")
    ]
    plot1.layer[4]["encoding"]["color"]["field"]="gen_Percent_Loading"
    plot1.layer[4]["encoding"]["color"]["scale"]["domain"]=[0,100]
    plot1.layer[4]["encoding"]["color"]["title"]="Gen Utilization %"
    plot1.layer[4]["encoding"]["size"]=Dict("field"=>"GenPower", "title"=>"Gen BaseMW", "type"=>"quantitative", "scale"=>Dict("range"=>[20,100]))
    

    # Layer 1 refers to branch data, layer 3 to bus data and layer 4 to generators data. 
    plot1.layer[1]["encoding"]["color"]["legend"]= Dict("orient"=>"bottom-right", "offset"=>-30)
    plot1.layer[1]["encoding"]["color"]["title"] = "Branch Utilization %"
    plot1.layer[1]["encoding"]["color"]["scale"]["domain"]= [0,100]

    if bus_data_input == "flex"
        plot1.layer[3]["encoding"]["color"]["legend"]= Dict("orient"=>"bottom-right")
        plot1.layer[3]["encoding"]["color"]["title"]= "Flexibility offered %"
        plot1.layer[3]["encoding"]["color"]["scale"]["domain"]= [0,20]

    elseif bus_data_input == "vm"
        plot1.layer[3]["encoding"]["color"]["legend"]= Dict("orient"=>"bottom-right")
        plot1.layer[3]["encoding"]["color"]["title"]= "Voltage magnitude"
        v_max, v_min = dict_find(net_data, "bus", "vm")
        plot1.layer[3]["encoding"]["color"]["scale"]["domain"]= [v_min,v_max]
    
    end
    

    plot1.layer[4]["encoding"]["color"]["legend"]=Dict("orient"=>"bottom-right")

    @set! plot1.resolve.scale.size=:independent
    @set! plot1.resolve.scale.color=:shared

    if zoom
        PowerPlots.Experimental.add_zoom!(plot1)
    end

    plot1


end

plot_grid(net_data, "vm", false)