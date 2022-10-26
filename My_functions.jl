#GENERATE VOLTAGE PROFILE FOR ALL FEEDERS FROM THE INPUT FILE
#RELEVANT DATA STORED IN feeder_ID

include("Network_paths.jl")
using ColorSchemes

#####################################################################################################################

## FUNCTIONS TO BE CALLED BEFORE PF ##

#####################################################################################################################


# Insert in net_data the flexibility % offered by each load, needed in my_build_pf
function add_load_flexibility(net_data::Dict, load_flexibility::Float64)

    flex = Dict{String,Any}()
    [flex[i] = Dict("flex_%" => load_flexibility) for (i,load) in net_data["load"]]

    return Dict{String,Any}("load" => flex)

end

# Insert in net_data the capacity congested (in %) of each branch, needed in my_build_pf
function add_congestion_capacity(net_data::Dict, cap::Float64)

    congestion_capacity = Dict{String,Any}()
    [congestion_capacity[i] = Dict("cong_cap" => cap) for (i,branch) in net_data["branch"]]

    return Dict{String,Any}("branch" => congestion_capacity)
end

# Add generator with curtailment capabilities 
function add_single_generator(net_data::Dict, size, gen, curt)

    if isa(gen, String)
        gen = parse(Int64, gen)
    end
    
    i = length(net_data["gen"]) + 1

    net_data["gen"]["$i"] = Dict("pg" =>size, "qg" =>0, "pmin" => size * (1-curt), "pmax"=>size, "qmin" =>0, "qmax"=>0, "gen_bus" => gen, "gen_status"=>1, "index" => i, "source_id" => ["gen", i])
    net_data["bus"]["$gen"]["bus_type"] = 2
    
end


#####################################################################################################################

## PF FUNCTIONS ##

#####################################################################################################################


# Add congestion constraints iteratively: reduce computational time 
function solve_pf_branch_power_cuts_mine(data::Dict{String,<:Any}, model_type::Type, optimizer, ref; kwargs...)
    return solve_pf_branch_power_cuts_mine!(data, model_type, optimizer, ref; kwargs...)
end

function solve_pf_branch_power_cuts_mine!(data::Dict{String,<:Any}, model_type::Type, optimizer, ref; solution_processors=[], max_iter::Int=100, time_limit::Float64=3600.0)

    Memento.info(_LOGGER, "maximum cut iterations set to value of $max_iter")

    for (i,branch) in data["branch"]
        if haskey(branch, "rate_a")
            branch["rate_a_inactive"] = branch["rate_a"]
        end
    end

    start_time = time()

    pm = instantiate_model(data, model_type, ref)
    result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=solution_processors)

    iteration = 1
    violated = true
    while violated && iteration < max_iter && (time() - start_time) < time_limit
        violated = false

        # Congestion constraint
        for (i,branch) in data["branch"]
            if haskey(branch, "rate_a_inactive")
                rate_a = branch["rate_a_inactive"]
                branch_sol = result["solution"]["branch"][i]

                cap = branch["cong_cap"]

                mva_fr = abs(branch_sol["pf"])
                mva_to = abs(branch_sol["pt"])

                if !isnan(branch_sol["qf"]) && !isnan(branch_sol["qt"])
                    mva_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                    mva_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)
                end

                #println(branch["index"], rate_a, mva_fr, mva_to)

                if mva_fr > rate_a * cap || mva_to > rate_a * cap
                    Memento.info(_LOGGER, "activate rate_a * cap on branch $(branch["index"])")

                    branch["rate_a"] = branch["rate_a_inactive"]
                    delete!(branch, "rate_a_inactive")

                    idx = branch["index"]
                    constraint_thermal_limit_from_me(pm, idx)
                    constraint_thermal_limit_to_me(pm, idx)

                    violated = true
                end
            end
        end

        if violated
            iteration += 1
            result = optimize_model!(pm, solution_processors=solution_processors)
        else
            Memento.info(_LOGGER, "flow cuts converged in $iteration iterations")
        end
    end

    result["solve_time"] = time() - start_time
    result["iterations"] = iteration

    return result
end

# Add congestion and voltage constraints ieratively when violated: reduce computational time 
function solve_pf_branch_voltage_cuts_mine(data::Dict{String,<:Any}, model_type::Type, optimizer,ref; kwargs...)
    return solve_pf_branch_voltage_cuts_mine!(data, model_type, optimizer, ref; kwargs...)
end

function solve_pf_branch_voltage_cuts_mine!(data::Dict{String,<:Any}, model_type::Type, optimizer, ref; solution_processors=[], max_iter::Int=100, time_limit::Float64=3600.0)

    Memento.info(_LOGGER, "maximum cut iterations set to value of $max_iter")

    for (i,branch) in data["branch"]
        if haskey(branch, "rate_a")
            branch["rate_a_inactive"] = branch["rate_a"]
        end
    end

    start_time = time()

   
    pm = instantiate_model(data, model_type, ref)
    result = optimize_model!(pm, optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0), solution_processors=solution_processors)

    #print_summary(result["solution"])

    iteration = 1
    violated = true
    while violated && iteration < max_iter && (time() - start_time) < time_limit
        violated = false

        # Congestion constraint
        for (i,branch) in data["branch"]
            if haskey(branch, "rate_a_inactive")
                rate_a = branch["rate_a_inactive"]
                branch_sol = result["solution"]["branch"][i]

                cap = branch["cong_cap"]

                mva_fr = abs(branch_sol["pf"])
                mva_to = abs(branch_sol["pt"])

                if !isnan(branch_sol["qf"]) && !isnan(branch_sol["qt"])
                    mva_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                    mva_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)
                end

                #println(branch["index"], rate_a, mva_fr, mva_to)

                if mva_fr > rate_a * cap || mva_to > rate_a * cap
                    Memento.info(_LOGGER, "activate rate_a * cap on branch $(branch["index"])")

                    branch["rate_a"] = branch["rate_a_inactive"]
                    delete!(branch, "rate_a_inactive")

                    idx = branch["index"]
                    constraint_thermal_limit_from_me(pm, idx)
                    constraint_thermal_limit_to_me(pm, idx)

                    violated = true
                end
            end
        end

        # Voltage constraint 
        for (i,bus) in data["bus"]

            vm = round(result["solution"]["bus"][i]["vm"],digits = 4)
            v_min = bus["vmin"]
            v_max = bus["vmax"]
            idx = bus["index"]
            if vm < v_min
    
                Memento.info(_LOGGER, "activate v_min on bus $i")
                constraint_voltage_magnitude_lower_me(pm,idx)
                
                violated = true
            end

            if vm > v_max

                Memento.info(_LOGGER, "activate v_max on bus $i")
                constraint_voltage_magnitude_upper_me(pm,idx)

                violated = true
                
            end
        end

        if violated
            iteration += 1
            #result = solve_opf(data, model_type, optimizer)
            result = optimize_model!(pm, solution_processors=solution_processors)

            #print_summary(result["solution"])
        else
            Memento.info(_LOGGER, "flow cuts converged in $iteration iterations")
        end
    end

    result["solve_time"] = time() - start_time
    result["iterations"] = iteration

    return result
end









#####################################################################################################################

## FUNCTIONS FOR NETWORK ANALYSIS ONCE PF IS DONE ##

#####################################################################################################################


# Compute & plot voltage profile of each feeder
function calc_voltage_profile(net_data::Dict, result::Dict, file_name, do_plot = false, save_fig = false, boundaries = false)
   
    # INPUT FOR PATH CALC

    down_node, g, map = downstreamcalcs(net_data)

    extremes = []
    [push!(extremes, i) for (i,j) in down_node if j ==[]]

    tot_paths = [Tuple[]]
    dist = AbstractFloat[]

    len = Dict((map[j["from_b"]],map[j["to_b"]]) => j["length"] for (ind,j) in net_data["distance"])

    inv_map = inverse_map(net_data)
    node_dist = []

    tot_paths, dist, node_dist = paths_calc(g, extremes, len, dist, tot_paths, inv_map, file_name, node_dist)

    paths = []
    feeders = []

    mv_busbar = 0

    #COUNT NUMBEER OF FEEDERS
    for (path,line) in enumerate(tot_paths)
    
        if line[1][2]!=1 && line[1][2] ∉ feeders
            push!(feeders, line[1][2])
        end
    end

    # CREATE DICT CONTAINING ALL DATA RELEVANT FOR FEEDERS
    alphabet = 'A':'Z'
    feeder_ID = Dict( j => Dict("Name" => "Feeder $i", "Paths" => [], "Paths_ID" => [], "Paths_distance"=> [],"Paths_volt" =>[],"vmin" =>[],"vmax" =>[]) for (i,j) in zip('A':alphabet[length(feeders)], feeders) )

    # Write tot_paths in a vectorial way. 
    for vec in tot_paths

        feed = []
        for j in 1:length(vec)

            push!(feed,vec[j][1])
            
        end
    
        push!(paths,feed)  #Rewrite how paths are made

    end

    for (idx,path) in enumerate(paths)

        if path[end]!=1 && path[2] in keys(feeder_ID)
            push!(feeder_ID[path[2]]["Paths"], path)
            push!(feeder_ID[path[2]]["Paths_ID"], idx)
            push!(feeder_ID[path[2]]["Paths_distance"], node_dist[idx])
        end

    end

    # Volt profile of each path

    path_volt = sort(Dict(i => zeros(length(tot_paths[i])) for i in keys(tot_paths)))

    # GET VECTOR OF VOLTAGE PROFILES

    for (path_ID, buses) in zip(keys(path_volt),paths)

        for (bus, volt) in result["solution"]["bus"]
            
            bus = parse(Int64,bus)
            if bus in buses
                ind = findall(x-> x==bus, buses)[1]
                path_volt[path_ID][ind] = volt["vm"]
            end
        end
    end

    #= ANALYSIS ON VOLT PROFILE
    lowest_voltage = Dict( i => Dict() for i in keys(path_volt) )

    min = [1.1]
    critical_line = [0]
    for (path, feed) in path_volt

        # Find node with lowest voltage 
        value, bus_ind = findmin(feed)
        lowest_voltage[path] = Dict(paths[path][bus_ind] => value) 
        println("For path $path the lowest voltage is at node $bus_ind with $value p.u.\n")

        if value < min[1] && path!=1   # does not make sense to consider path 1 

            min[1] = value
            critical_line[1] = path

        end
    end
    =#
    for feeder in feeders
        [push!(feeder_ID[feeder]["Paths_volt"],path_volt[i]) for (j,i) in enumerate(feeder_ID[feeder]["Paths_ID"])]
    end
    
    # Add ID of branches belonging to the same feeder

    for (ref, data) in feeder_ID

        branches = []
    
        for path in feeder_ID[ref]["Paths"]
            
            for (ind,j) in enumerate(path)
        
                if ind<length(path)
        
                    t_bus = path[ind]
                    f_bus = path[ind+1]
        
                    for (idx,branch) in net_data["branch"]
        
                        if f_bus == branch["f_bus"] && t_bus == branch["t_bus"]
                            if idx ∉ branches
                                push!(branches,idx)
                            end
                        end
                    end
                end
            end
        end
    
        feeder_ID[ref]["Branches"] = branches
    
    end

    # Add ID of bus belonging to the same feeder
    for (ref, data) in feeder_ID

        nodes = []

        for path in feeder_ID[ref]["Paths"]
            [push!(nodes, bus) for bus in path if bus ∉ nodes]
        end

        feeder_ID[ref]["Buses"] = nodes
        mv_busbar = nodes[1]
    end

    for feeder in feeders

        nrows = maximum(extrema(length, values(feeder_ID[feeder]["Paths_distance"])))
        ncols = length(feeder_ID[feeder]["Paths_distance"])

        y = fill(NaN, nrows, ncols)
        x = fill(NaN, nrows, ncols)

        [y[1:length(path_volt[i]),j] = path_volt[i] for (j,i) in enumerate(feeder_ID[feeder]["Paths_ID"])]
        [x[1:length(i),j] = i for (j,i) in enumerate(feeder_ID[feeder]["Paths_distance"])]
        
        min = minimum(filter(!isnan,y))
        max = maximum(filter(!isnan,y))

        feeder_ID[feeder]["vmin"] = min
        feeder_ID[feeder]["vmax"] = max

        feeder_ID

        plot = Plots.plot(x,y; marker=(:circle,4), linewidth = 2, label = "")
        Plots.plot!([min], color=:black ,seriestype = "hline", linewidth = 3, label = "v_min = $min pu", legend = :bottomleft, left_margin = 5Plots.mm, right_margin = 15Plots.mm)
        Plots.plot!([max], color=:black ,seriestype = "hline", linewidth = 3, label = "v_max = $max pu", legend = :bottomleft)

        if boundaries
            Plots.plot!([0.9], color=:red ,seriestype = "hline", linewidth = 3, label = "v = 0.9 pu")
            Plots.plot!([1.1], color=:red ,seriestype = "hline", linewidth = 3, label = "v = 1.1 pu")
        end
        
        xlabel!("Distance (km)")
        ylabel!("Voltage drop (p.u.)")
        title!("Voltage drop for: $(feeder_ID[feeder]["Name"])")
        Plots.plot!(size = (1000,800))

        if do_plot
            display(plot)
        end

        if save_fig

            file = replace(file_name, "Official_"=>"")
            file = replace(file,"_test.m" =>"")
            file = uppercasefirst(file)
            
            Plots.savefig("C:\\Users\\u0152683\\Desktop\\Networks\\PF simulation\\Flexible nodes\\DGs\\Rural\\Voltage profile\\P=14\\$(feeder_ID[feeder]["Name"]).png")
        end
    end
    

    return feeder_ID, path_volt, mv_busbar

end

# Computes loading of each branch (adding it to net_data) and updates gen_ID
function calc_branch_loading(net_data::Dict, feeder_ID::Dict,gen_ID::Dict, congestion_limit)
   
    [net_data["branch"][i]["loading"] = round(max(abs(complex(branch["pt"],branch["qt"]))/branch["rate_a"],abs(complex(branch["pf"],branch["qf"]))/branch["rate_a"]), digits = 3)*100 for (i,branch) in net_data["branch"]]
    
    congested_lines = []
    for (i,data) in net_data["branch"]
        if data["loading"] > congestion_limit
           push!(congested_lines,i)  
        end
    end

    for (idx,gen) in gen_ID

        branch_load = []
        ref = gen["ref"]
    
        [push!(branch_load,net_data["branch"][i]["loading"]) for i in feeder_ID[ref]["Branches"]]
    
        max_load, ind_max = findmax(branch_load)
        gen_ID[idx]["max_branch_load"] = max_load
        gen_ID[idx]["Critical_branch"] = feeder_ID[ref]["Branches"][ind_max]
        
    end

    for (id, data) in feeder_ID
        loadings = []
        [push!(loadings, net_data["branch"][i]["loading"]) for i in feeder_ID[id]["Branches"]]
        feeder_ID[id]["Branch_loading"] = loadings

        if any(in(data["Branches"]).(congested_lines)) 
            feeder_ID[id]["Congested_lines"] = congested_lines[in(data["Branches"]).(congested_lines)]
        else
            feeder_ID[id]["Congested_lines"] = []
        end
    end
end

# Computes power losses for each branch. 
function calc_power_losses(data::Dict{String,<:Any})
    
    pm_data = get_pm_data(data)

    losses = Dict{String,Any}()
    p_loss = []
    q_loss = []

    for (i,branch) in pm_data["branch"]

        if branch["br_status"] != 0
            pt = branch["pt"]
            qt = branch["qt"]
            pf = branch["pf"]
            qf = branch["qf"]

            losses[i] = Dict(
                "p_loss" => pt+pf,
                "q_loss" => qt+qf
            )
            push!(p_loss, pt+pf)
            push!(q_loss, qt+qf)
        else
            losses[i] = Dict(
                "p_loss" => NaN,
                "q_loss" => NaN
            )
            push!(p_loss, pt+pf)
            push!(q_loss, qt+qf)
        end

    end

    return Dict{String,Any}("branch" => losses), p_loss, q_loss

end

# Evaluate the flexibility offered by each load
function calc_flexibility_offered(net_data::Dict, result::Dict)

    flex_loads = Dict()
    p_load = 0
    q_load = 0

    for (i,load) in result["solution"]["load"]
        p_load += load["load_p"]
        q_load += load["load_q"]
        

        if load["load_p"] < net_data["load"][i]["pd"] || load["load_q"] < net_data["load"][i]["qd"]

            diff_p = net_data["load"][i]["pd"]-load["load_p"]
            diff_q = net_data["load"][i]["qd"]-load["load_q"]
            
            flex_loads[string(net_data["load"][i]["load_bus"])] = Dict(
                "diff_real" => diff_p,
                "diff_imm" => diff_q,
                "flex_p" => round((diff_p)/net_data["load"][i]["pd"], digits = 3)*100,
                "flex_q" => round((diff_q)/net_data["load"][i]["qd"], digits = 3)*100,
            )
            net_data["bus"][string(net_data["load"][i]["load_bus"])]["flex"] = round((diff_p)/net_data["load"][i]["pd"], digits = 3)*100

        else
            
            flex_loads[string(net_data["load"][i]["load_bus"])] = Dict(
                "diff_real" => 0,
                "diff_imm" => 0,
                "flex_p" => 0,
                "flex_q" => 0,
            )
            net_data["bus"][string(net_data["load"][i]["load_bus"])]["flex"] = 0
        end

    end

    return flex_loads, p_load, q_load

end

# Evaluate flexibility offered by each load taking into account Upwards and Downwards flexbility 
function calc_flexibility_offered_new(net_data::Dict, result::Dict)
    flex_loads = Dict()
    p_load = sum([load["load_p"] for (x,load) in result["solution"]["load"]])
    q_load = sum([load["load_q"] for (x,load) in result["solution"]["load"]])

    for (i,load) in result["solution"]["load"]

        
        diff_p = load["load_p"] - net_data["load"][i]["pd"]
        diff_q = load["load_q"] - net_data["load"][i]["qd"]

        if diff_p < 0 && abs(diff_p) > 10^-5

            flex_loads[string(net_data["load"][i]["load_bus"])] = Dict(
                "diff_real" => diff_p,
                "flex_p" => round((diff_p)/net_data["load"][i]["pd"], digits = 3)*100,
            )
            net_data["bus"][string(net_data["load"][i]["load_bus"])]["p_flex"] = round((diff_p)/net_data["load"][i]["pd"], digits = 3)*100
            net_data["bus"][string(net_data["load"][i]["load_bus"])]["flex_type"] = "Downward"

        elseif diff_p > 0 && abs(diff_p) > 10^-5

            flex_loads[string(net_data["load"][i]["load_bus"])] = Dict(
                "diff_real" => diff_p,
                "flex_p" => round((diff_p)/net_data["load"][i]["pd"], digits = 3)*100,
            )
            net_data["bus"][string(net_data["load"][i]["load_bus"])]["p_flex"] = round((diff_p)/net_data["load"][i]["pd"], digits = 3)*100
            net_data["bus"][string(net_data["load"][i]["load_bus"])]["flex_type"] = "Upward"
        else
            flex_loads[string(net_data["load"][i]["load_bus"])] = Dict(
                "diff_real" => 0.0,
                "flex_p" => 0,
            )
            net_data["bus"][string(net_data["load"][i]["load_bus"])]["p_flex"] = 0
            net_data["bus"][string(net_data["load"][i]["load_bus"])]["flex_type"] = "None"
        end
        
        if diff_q < 0 && abs(diff_q) > 10^-5

            flex_loads[string(net_data["load"][i]["load_bus"])]["diff_imm"] = diff_q
            flex_loads[string(net_data["load"][i]["load_bus"])]["flex_q"] = round((diff_q)/net_data["load"][i]["qd"], digits = 3)*100

            net_data["bus"][string(net_data["load"][i]["load_bus"])]["q_flex"] = round((diff_q)/net_data["load"][i]["qd"], digits = 3)*100
        
        elseif diff_q > 0 && abs(diff_q) > 10^-5

            flex_loads[string(net_data["load"][i]["load_bus"])]["diff_imm"] = diff_q
            flex_loads[string(net_data["load"][i]["load_bus"])]["flex_q"] = round((diff_q)/net_data["load"][i]["qd"], digits = 3)*100

            net_data["bus"][string(net_data["load"][i]["load_bus"])]["q_flex"] = round((diff_q)/net_data["load"][i]["qd"], digits = 3)*100
        else
            flex_loads[string(net_data["load"][i]["load_bus"])]["diff_imm"] = 0.0
            flex_loads[string(net_data["load"][i]["load_bus"])]["flex_q"] = 0.0
            
            net_data["bus"][string(net_data["load"][i]["load_bus"])]["q_flex"] = 0
        end
        
        v = []
        [push!(v,string(data["load_bus"])) for (ID,data) in net_data["load"]]
        b = collect(keys(net_data["bus"]))
        miss = findall(x->x==0, in(v).(b))
        if !isempty(b[miss])
            for idx in b[miss]

                flex_loads[idx] = Dict(
                    "diff_real" => 0.0,
                    "diff_imm" => 0.0,
                    "flex_p" => 0,
                    "flex_q" => 0,
                )
                net_data["bus"][idx]["p_flex"] = 0
                net_data["bus"][idx]["q_flex"] = 0
            end
        end

    end



    return flex_loads, p_load, q_load

end

# Evaluate flexibility offered by each load taking into account Upwards and Downwards flexbility 
function calc_flexibility_offered_final(net_data::Dict, result::Dict)

    flex_loads = Dict()

    p_load = sum([load["load_p"] for (x,load) in result["solution"]["load"]])
    q_load = sum([load["load_q"] for (x,load) in result["solution"]["load"]])


    for (i,load) in result["solution"]["load"]

        p_nominal = net_data["load"][i]["pd"]
        q_nominal = net_data["load"][i]["qd"]

        bus_id = string(net_data["load"][i]["load_bus"])
        flex_loads[bus_id] = Dict("diff_real" => 0.0, "flex_p" => 0.0, "diff_imm" => 0.0, "flex_q" => 0.0)

        net_data["bus"][bus_id]["p_flex"] = 0
        net_data["bus"][bus_id]["q_flex"] = 0
        net_data["bus"][bus_id]["flex_type"] = "None"
        
        x_p = load["x_p"]
        y_p = load["y_p"]

        x_q = load["x_q"]
        y_q = load["y_q"]

        if x_p > 10^-6  # upwards active flexibiility 

            flex_loads[bus_id]["diff_real"] = x_p
            flex_loads[bus_id]["flex_p"] = round(x_p/p_nominal, digits = 3)*100

            net_data["bus"][bus_id]["p_flex"] = flex_loads[bus_id]["flex_p"]
            net_data["bus"][bus_id]["flex_type"] = "Upward"
        
        elseif y_p > 10^-6

            flex_loads[bus_id]["diff_real"] = - y_p
            flex_loads[bus_id]["flex_p"] = - round(y_p/p_nominal, digits = 3)*100

            net_data["bus"][bus_id]["p_flex"] = flex_loads[bus_id]["flex_p"]
            net_data["bus"][bus_id]["flex_type"] = "Downward"
        end

        if x_q > 10^-6  # upwards reactive flexibiility 

            flex_loads[bus_id]["diff_imm"] = x_q
            flex_loads[bus_id]["flex_q"] = round(x_q/q_nominal, digits = 3)*100

            net_data["bus"][bus_id]["q_flex"] = flex_loads[bus_id]["flex_q"]
            net_data["bus"][bus_id]["flex_type"] = "Upward"
        
        elseif y_q > 10^-6

            flex_loads[bus_id]["diff_imm"] = - y_q
            flex_loads[bus_id]["flex_q"] = - round(y_q/q_nominal, digits = 3)*100

            net_data["bus"][bus_id]["q_flex"] = flex_loads[bus_id]["flex_q"]
            net_data["bus"][bus_id]["flex_type"] = "Downward"
        end

        
        v = []
        [push!(v,string(data["load_bus"])) for (ID,data) in net_data["load"]]
        b = collect(keys(net_data["bus"]))
        miss = findall(x->x==0, in(v).(b))
        if !isempty(b[miss])
            for idx in b[miss]

                flex_loads[idx] = Dict(
                    "diff_real" => 0.0,
                    "diff_imm" => 0.0,
                    "flex_p" => 0,
                    "flex_q" => 0,
                )
                net_data["bus"][idx]["p_flex"] = 0
                net_data["bus"][idx]["q_flex"] = 0
            end
        end

    end



    return flex_loads, p_load, q_load

end










#####################################################################################################################

## FUNCTIONS OF GENERAL INTEREST ##

#####################################################################################################################


# relating gen position to feeder
function get_gen_info(net_data::Dict, feeder_ID::Dict)

    generator_ID = Dict()

    for (j, feeder) in feeder_ID

        for (i,gen) in net_data["gen"]
            
            present = false

            for k in feeder["Paths"]
                if gen["gen_bus"] in k
                    present = true
                end
            end

            if present && gen["gen_bus"]!="1"  #exclude slack generator
                generator_ID[i] = Dict("ref" => j, "feeder"=> feeder["Name"], "bus"=>gen["gen_bus"])
            end
        end
    end


    return generator_ID

end

# Find maximum or min of certain proprieties in a dict (i.e. max voltage etc) 
function dict_find(d::Dict, category::String, property::String, max::Bool=true, min::Bool=true)

    vector = []
    try
        [push!(vector,data[property]) for (ID, data) in d[category]]
    catch e
        if isa(e, KeyError)
            println("Invalid category or property: ", e)
            println("Re-insert category:")
            cat = readline()
            println("Re-insert property: ")
            prop = readline()
           
            return dict_find(d, cat, prop )

        end
    end

    if max && min
        return maximum(vector), minimum(vector)

    elseif max
        return maximum(vector), nothing

    elseif min
        return nothing, minimum(vector)
        
    else
        println("Error: Select at least a max or as min ")
        return nothing, nothing 
    end

end

function dict_find_new(d::Dict, category::String, property::String, max::Bool=true, min::Bool=true)

    vector = []
    
    [push!(vector,get(data,property,0)) for (ID, data) in d[category]]

    if max && min
        return maximum(vector), minimum(vector)

    elseif max
        return maximum(vector), nothing

    elseif min
        return nothing, minimum(vector)
        
    end

end

# Find the node/s with maximum PTDF for a certain branch 
function calc_branch_max_ptdf(net_data::Dict, branch)

    if typeof(branch) != Int64
        branch = parse(Int64, branch)
    end

    # Convert from actual value to basic structure for ptdf computation
    basic_data = make_basic_network(deepcopy(net_data))
    ptdf = calc_basic_ptdf_matrix(basic_data)

    max = maximum(ptdf[branch,:])
    idx = [i for (i, x) in enumerate(ptdf[branch,:]) if x == max]
    return max, [basic_to_real_bus(net_data)[i] for i in idx]
end

# Get the value of the PTDF of a specific node on a specific branch 
function calc_node_branch_ptdf(net_data::Dict, branch, node)

    if typeof(branch) != Int64
        branch = parse(Int64, branch)
    end

    if typeof(node) != Int64
        node = parse(Int64, node)
    end

    node = real_to_basic_bus(net_data)[node]

    # Convert from actual value to basic structure for ptdf computation
    basic_data = make_basic_network(deepcopy(net_data))
    ptdf = calc_basic_ptdf_matrix(basic_data)

    return ptdf[branch,node]
end

# get the ptdf of a bus ( you can pass a vector of busses) on a certain branch in terms or reactive or active power. 
# ptdf computed by solving PF before and after the power variation  
function calc_my_ptdf(net_data::Dict, bus::Vector, br::Int64, var::Float64, type = "active")

    function sensitivity(net_data::Dict, ref, data::Dict, flows, bus, br::Int64, var::Float64, type)

        bus = ref[:bus_loads][bus][1] 

        if type== "active"

            data["load"]["$bus"]["pd"] = net_data["load"]["$bus"]["pd"]+var
        
            result_2 = solve_ac_pf(data, Ipopt.Optimizer)
            update_data!(data, result_2["solution"])
            flows_2 = calc_branch_flow_ac(data)
            update_data!(data, flows_2)
        
            ptdf_from = (flows_2["branch"]["$br"]["pf"] - flows["branch"]["$br"]["pf"])/var
            ptdf_to = (flows_2["branch"]["$br"]["pt"] - flows["branch"]["$br"]["pt"])/var

            return ptdf_from, ptdf_to 
        elseif type =="reactive"

            data["load"]["$bus"]["qd"] = net_data["load"]["$bus"]["qd"]+var
        
            result_2 = solve_ac_pf(data, Ipopt.Optimizer)
            update_data!(data, result_2["solution"])
            flows_2 = calc_branch_flow_ac(data)
            update_data!(data, flows_2)
        
            ptdf_from = (flows_2["branch"]["$br"]["qf"] - flows["branch"]["$br"]["qf"])/var
            ptdf_to = (flows_2["branch"]["$br"]["qt"] - flows["branch"]["$br"]["qt"])/var

            return ptdf_from, ptdf_to
        end

    
    end

    ref = PowerModels.build_ref(net_data)[:it][:pm][:nw][0]

    result = solve_ac_pf(net_data, Ipopt.Optimizer)
    update_data!(net_data, result["solution"])
    flows = calc_branch_flow_ac(net_data)
    update_data!(net_data, flows)

    ptdf = Dict()

    for b in bus

        data = deepcopy(net_data)

        ptdf_from, ptdf_to = sensitivity(net_data,ref,data, flows, b,br,var, type)  #6920, 81, -0.01
        ptdf[b] = Dict("ptdf_from" => ptdf_from,"ptdf_to" => ptdf_to, "type" => type)

    end

    return ptdf
end

# print different things
function printing_statements(result, p_load, q_load, p_loss, q_loss, flex_loads, gen_ID, flows, index, mv_busbar)

    println("\n Termination status: ", result["termination_status"])
    print("\n Objective function: ", result["objective"])

    print("\n Real power generated: ", result["solution"]["gen"]["1"]["pg"], "MW \n")
    print("\n Immaginary power generated: ", result["solution"]["gen"]["1"]["qg"], "MVar \n")

    print("\n MV busbar voltage : ", round(result["solution"]["bus"]["$mv_busbar"]["vm"],digits =5), "p.u. \n")

    print("\n Real power requested: ", p_load, "MW \n")
    print("\n Immaginary power requested: ", q_load, "MVar \n")

    print("\n Real power losses: ", sum(p_loss), "MW \n")
    print("\n Immaginary power losses: ", sum(q_loss), "MVar \n")

    print("\n Active power Flexibility offered: ", sum([abs(flex["diff_real"]) for (i,flex) in flex_loads]), "MW \n")
    print("\n Rective power Flexibility offered: ", sum([abs(flex["diff_imm"]) for (i,flex) in flex_loads]), "MVar \n")

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

    #[println("\n Vmin: ", round(feeder_ID[gen["ref"]]["vmin"], digits = 5)," Vmax: ",round(feeder_ID[gen["ref"]]["vmax"],digits = 5)," for ",gen["feeder"]," (where generator $i is connected)") for (i,gen) in gen_ID]
    #[println("\n Feeder loading: ", round(gen["max_branch_load"],digits = 3),"% in branch ",gen["Critical_branch"]) for (i,gen) in gen_ID]

    feeder_names = Vector{String}()
    [push!(feeder_names, feeder["Name"]) for (i,feeder) in feeder_ID]
    sort!(feeder_names)

    for (id, gen) in gen_ID
        println(" Generator $id is connected at bus ", gen["bus"], " in ", gen["feeder"])
    end

    println()

    for f_name in feeder_names

        for (ref,data) in feeder_ID
            if data["Name"] == f_name

                println(" ",f_name," --> Vmin: ", round(data["vmin"], digits = 5), " Vmax: ", round(data["vmin"], digits = 5))
                
            end
        end
    end
    println()



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
    
    num = length(filter(x->x>index, branch_loading))

    println("\nBranches with loading higher than ", index, "%: ", num)


    vmin = Vector()
    vmax = Vector()
    [push!(vmin,feeder["vmin"]) for (id, feeder) in feeder_ID ]
    println("\n Lowest voltage magnitude in the grid: ", minimum(vmin))
    [push!(vmin,feeder["vmax"]) for (id, feeder) in feeder_ID ]
    println("\n Highest voltage magnitude in the grid: ", maximum(vmin))

end

# Create Dataframe to store data on feeders concerning congestion issues and voltage issues. 
# If save=true then the data is stored in the xlsx file 
function create_df(feeder_ID, result, threshold, save = false)
    
    df = DataFrame()

    feeder_names = Vector{String}()
    [push!(feeder_names, feeder["Name"]) for (i,feeder) in feeder_ID]
    sort!(feeder_names)

    current_violation = Vector()
    c_v_value = Vector()
    c_v_branch = Vector()
    c_v_max = Vector()

    voltage_violation = Vector()
    v_v_bus = Vector()
    v_v_type = Vector()
    v_min = Vector()
    v_max = Vector()

    for f_name in feeder_names
        
        upper = false
        lower = false

        for (id,data) in feeder_ID
            if data["Name"] == f_name
                # Congestion analysis 
                if !isempty(data["Congested_lines"])
                    push!(current_violation,true) 
                    push!(c_v_value, join(filter(x->x>threshold, data["Branch_loading"]), " "))
                    push!(c_v_branch, join(data["Congested_lines"], " "))
                else
                    push!(current_violation,false)
                    push!(c_v_value, "-")
                    push!(c_v_branch, "-")
                end
                
                v = Vector()
                #Voltage analysis
                for (i,bus) in result["solution"]["bus"]
                    i = parse(Int64,i)
                    if i in data["Buses"]
                        if bus["vm"]>1.1
                            push!(v, i)
                            upper = true
                        elseif bus["vm"]<0.9
                            push!(v, i)
                            lower = true
                        end
                    end
                end

                if !isempty(v)
                    push!(v_v_bus, join(v, " "))
                else
                    push!(v_v_bus, "-")
                end
                
                if upper && lower
                    push!(v_v_type,"UP, DWN")
                    push!(voltage_violation,true)
                elseif upper
                    push!(v_v_type,"UP")
                    push!(voltage_violation,true)
                elseif lower
                    push!(v_v_type,"DWN")
                    push!(voltage_violation,true)
                else
                    push!(v_v_type,"-")
                    push!(voltage_violation,false)
                end

                push!(v_min, data["vmin"])
                push!(v_max, data["vmax"])
                push!(c_v_max, maximum(data["Branch_loading"]))
            end

        end
    end

    df.feeders = replace.(feeder_names,"Feeder "=>"")
    df.Current_violation = current_violation
    df.Overloaded_branches = c_v_branch
    df.Overload = c_v_value
    df.Max_loading = c_v_max
    df.Voltage_violation = voltage_violation
    df.Violation_bus = v_v_bus
    df.Violation_type = v_v_type
    df.Vmin = round.(v_min,digits = 5)
    df.Vmax = round.(v_max, digits = 5)

    if save
        XLSX.openxlsx("C:\\Users\\u0152683\\Desktop\\Networks\\Comparison_networks.xlsx", mode="rw") do xf
            sheet = xf["Data"]
            for r in 1:size(df,1), c in 1:size(df,2)
                sheet[XLSX.CellRef(r,c)] = df[r,c]
            end
        end
    end

    return df

end











#####################################################################################################################

## PLOTTING FUNCTIONS ##

#####################################################################################################################




#Plot network with branch loading
# node_attribute can be either: "p_flex",  "q_flex", "vm", "flex_type","basic" 
# gen_attribute can be either: "pg", "curtailment", "qg", "basic"
# branch_attribute can be either: "loading", "q_loss", "p_loss", "pt", "basic"

function plot_grid(case::Dict, node_attribute::String, gen_attribute::String, branch_attribute::String; zoom::Bool=false, display_flow::Bool = true, save::Bool = false)

    # copy data for modification by plots
    data = deepcopy(case)

    # what components to display
    show_components = ["bus", "branch", "gen"]

    prop_node, prop_gen, prop_br = dict_of_proprieties(node_attribute, gen_attribute, branch_attribute)

    plot = powerplot(data, 
                    show_flow = display_flow,
                    components = show_components,
                    width = 1000, 
                    height = 1000,

                    branch_data = prop_br[:branch_data],  #branch
                    branch_data_type = prop_br[:branch_data_type],
                    branch_size = prop_br[:branch_size],
                    branch_color = prop_br[:color_range],

                    bus_data = prop_node[:bus_data],  #bus
                    bus_data_type = prop_node[:bus_data_type],
                    bus_size = prop_node[:bus_size],
                    bus_color = prop_node[:color_range],

                    gen_data = prop_gen[:gen_data],  #gen
                    gen_data_type = prop_gen[:gen_data_type],
                    gen_size = prop_gen[:gen_size],
                    gen_color = prop_gen[:color_range]
    )

    # BRANCH
    if branch_attribute != "basic"
        #plot.layer[1]["layer"][1]["encoding"]["color"]["field"]="branch_Percent_Loading"
        plot.layer[1]["layer"][1]["encoding"]["color"]["legend"]= Dict("orient"=>"bottom-right")
        plot.layer[1]["layer"][1]["encoding"]["color"]["title"]= prop_br[:title]
        plot.layer[1]["layer"][1]["encoding"]["color"]["scale"]["domain"]= prop_br[:range]
        plot.layer[1]["layer"][1]["encoding"]["color"]["scale"]["range"] = prop_br[:color_range]
        #plot.layer[1]["layer"][1]["encoding"]["size"]=Dict("field"=>"BranchPower", "title"=>"Branch BaseMW", "type"=>"quantitative", "scale"=>Dict("range"=>[3,10]))
    
    end

    # BUS
    if !(node_attribute in ["basic", "flex_type"])
        plot.layer[3]["encoding"]["color"]["legend"]= Dict("orient"=>"bottom-right", "offset" => -30)
        plot.layer[3]["encoding"]["color"]["title"]= prop_node[:title]
        plot.layer[3]["encoding"]["color"]["scale"]["domain"]= prop_node[:range]
        plot.layer[3]["encoding"]["color"]["scale"]["range"] = prop_node[:color_range]
    end
        
    # GEN
    if gen_attribute != "basic"
        #plot.layer[4]["encoding"]["color"]["field"]="gen_Percent_Loading"
        plot.layer[4]["encoding"]["color"]["legend"] = Dict("orient"=>"bottom-right", "offset" => -60)
        plot.layer[4]["encoding"]["color"]["title"] = prop_gen[:title]
        plot.layer[4]["encoding"]["color"]["scale"]["domain"] = prop_gen[:range]
        plot.layer[4]["encoding"]["color"]["scale"]["range"] = prop_gen[:color_range]
        #plot.layer[4]["encoding"]["size"]=Dict("field"=>"GenPower", "title"=>"Gen BaseMW", "type"=>"quantitative", "scale"=>Dict("range"=>[300,1000]))
    end
        
    @set! plot.resolve.scale.size=:independent
    #@set! plot.resolve.scale.color=:shared

    if zoom
        PowerPlots.Experimental.add_zoom!(plot)
    end

    plot


end

function dict_of_proprieties(node::String, gen::String, branch::String)
    
    dict_node = Dict{String, Dict}()
    dict_gen = Dict{String, Dict}()
    dict_branch = Dict{String, Dict}()

    f_p_max, f_p_min = dict_find_new(net_data, "bus", "p_flex")
    f_q_max, f_q_min = dict_find_new(net_data, "bus", "q_flex")
    v_max, v_min = dict_find_new(net_data, "bus", "vm")
    curt_max, min = dict_find_new(net_data, "gen", "curtailment")
    pg_max, pg_min = dict_find_new(net_data, "gen", "pg")
    qg_max, qg_min = dict_find_new(net_data, "gen", "qg")
    b_l_max, b_l_min = dict_find_new(net_data, "branch", "loading")
    p_loss_max, min = dict_find_new(net_data, "branch", "p_loss")
    q_loss_max, min= dict_find_new(net_data, "branch", "q_loss")
    p_f_max = dict_find_new(net_data, "branch", "pf")


    # Attributes for nodes (bus)
    dict_node["p_flex"] = Dict{Symbol, Any}(
        :bus_size => 90,
        :bus_data =>"p_flex",
        :bus_data_type => "quantitative",
        :range => [f_p_min, f_p_max],
        :title => "Active DR offered (%)",
        :color_range => ["#C0C0C0","#000000"]#colorscheme2array(ColorSchemes.colorschemes[:RdYlBu_10]),  #red - blue
    )
    dict_node["q_flex"] = Dict{Symbol, Any}(
        :bus_size => 90,
        :bus_data => "q_flex",
        :bus_data_type => "quantitative",
        :range => [f_q_min,f_q_max],
        :title => "Reactive DR offered (%)",
        :color_range => ["#C0C0C0","#000000"]#colorscheme2array(ColorSchemes.colorschemes[:PiYG_10]),  #pink - green
    )
    dict_node["vm"] = Dict{Symbol, Any}(
        :bus_size => 90,
        :bus_data => "vm",
        :bus_data_type => "quantitative",
        #:range => [v_min,v_max],
        :range => [0.9,1.1],
        :title => "Voltage Magnitude (p.u.)",
        :color_range => ["#FFF5EE","#B0E0E6","#000080"],
    )
    dict_node["flex_type"] = Dict{Symbol, Any}(
        :bus_size => 90,
        :bus_data => "flex_type",
        :bus_data_type => "nominal",
    )
    dict_node["basic"] = Dict{Symbol, Any}(
        :bus_size => 90,
        :bus_data => "ComponentType",
        :bus_data_type => "nominal",
        :color_range => "#228b22"  #ForestGreen
    )

    # Attributes for gen
    dict_gen["pg"] = Dict{Symbol, Any}(
        :gen_size => 150,
        :gen_data => "pg",
        :gen_data_type => "quantitative",
        :range => [pg_min,pg_max],
        :title => "Active power generated (MW)",
        :color_range => ["white","purple"],
    )
    dict_gen["qg"] = Dict{Symbol, Any}(
        :gen_size => 150,
        :gen_data => "qg",
        :gen_data_type => "quantitative",
        :range => [qg_min,qg_max],
        :title => "Reactive power generated (MVar)",
        :color_range => ["white","orange"],
    )
    dict_gen["curtailment"] = Dict{Symbol, Any}(
        :gen_size => 150,
        :gen_data => "curtailment",
        :gen_data_type => "quantitative",
        :range => [0,curt_max],
        :title => "DG Curtailment (%)",
        :color_range => colorscheme2array(ColorSchemes.colorschemes[:BuPu_3]),  #white to purple  
    )
    dict_gen["basic"] = Dict{Symbol, Any}(
        :gen_size => 150,
        :gen_data => "ComponentType",
        :gen_data_type => "nominal",
        :color_range => "purple"
    )

    # Attributes for branch
    dict_branch["loading"] = Dict{Symbol, Any}(
        :branch_size => 3,
        :branch_data => "loading",
        :branch_data_type => "quantitative",
        :range => [0,b_l_max],
        :title => "Branch loading (%)",
        :color_range => ["green","orange","red"],
    )
    dict_branch["q_loss"] = Dict{Symbol, Any}(
        :branch_size => 3,
        :branch_data => "q_loss",
        :branch_data_type => "quantitative",
        :range => [0,q_loss_max],
        :title => "Q loss (MVar)",
        :color_range => ["grey","red"],
    )
    dict_branch["p_loss"] = Dict{Symbol, Any}(
        :branch_size => 3,
        :branch_data => "p_loss",
        :branch_data_type => "quantitative",
        :range => [0,p_loss_max],
        :title => "P loss (MW)",
        :color_range => ["grey","red"],
    )
    dict_branch["pf"] = Dict{Symbol, Any}(
        :branch_size => 3,
        :branch_data => "pf",
        :branch_data_type => "quantitative",
        :range => [0,p_f_max],
        :title => "Power Flow (MW)",
        :color_range => ["green","red"],
    )
    dict_branch["basic"] = Dict{Symbol, Any}(
        :branch_size => 3,
        :branch_data => "ComponentType",
        :branch_data_type => "nominal",
        :color_range => "#87cefa" # lightskyblue #00BFF" #deepSkyBlue
    )

    return dict_node[node],dict_gen[gen], dict_branch[branch]

end

#Plot network with colored feeders
function plot_corall_grid(net_data, file_name, path_to_save,zoom = true, save = false  )

    feeder_ID_1, mv_busbar = get_feeder_data(net_data, file_name)
    title  = replace(file_name, "Official_"=>"")
    title = uppercasefirst(replace(title,".m" =>""))

    for (k, feeder) in feeder_ID_1
        for (id, branch) in net_data["branch"]
            if id in feeder["Branches"]
                net_data["branch"][id]["feeder"] = feeder["Name"]
            end
        end
        for (idx, bus) in net_data["bus"]
            if parse(Int64,idx) in feeder["Buses"][2:end]
                net_data["bus"][idx]["feeder"] = feeder["Name"]
            end
        end
    
    end
    
    [net_data["branch"][id]["feeder"] = "HV/MV" for (id,branch) in net_data["branch"] if !haskey(branch,"feeder")]
    [net_data["bus"][id]["feeder"] = "HV/MV" for (id,bus) in net_data["bus"] if !haskey(bus,"feeder")]
    
    if title != "Semiurban"
        color_range = colorscheme2array(ColorSchemes.colorschemes[:tableau_10])
        color_range = [color_range[i] for i in 1:length(feeder_ID_1)+2]
    else
        color_range = colorscheme2array(ColorSchemes.colorschemes[:tableau_20])
        color_range = [color_range[i] for i in 1:length(feeder_ID_1)+2]
    end
    
    plot = powerplot(net_data,width = 800, height = 800, branch_size = 4, bus_size = 15, gen_size = 100,
    branch_data = "feeder", branch_color = color_range,
    bus_data = "feeder", node_color = color_range,
    gen_color = "black" )
    
    
    @set! plot.title = title
    @set! plot.resolve.scale.color=:shared
    
    plot.layer[1]["encoding"]["color"]["title"]="Feeders"
    plot.layer[1]["encoding"]["color"]["legend"]=Dict("clipHeight"=>50, "type" => "symbol", "labelFontSize"=>10, "symbolType" => "circle", "symbolSize" => 1000)
    if zoom 
        PowerPlots.Experimental.add_zoom!(plot)
    end

    plot
end

# Update net_data with curtailment values for each gen
function calc_curtailment(net_data, result)

    for (i, gen) in result["solution"]["gen"]

        if i!="1"
            p_nominal = net_data["gen"][i]["pmax"]
            p_final = gen["pg"]
            curtail = p_nominal - p_final
        
            if curtail > 10^-4
                net_data["gen"][i]["curtailment"] = round(curtail/p_nominal, digits = 4)*100 
            else
                net_data["gen"][i]["curtailment"] = 0
            end
            
        else
            net_data["gen"][i]["curtailment"] = 0
        end
    end
    
end



#####################################################################################################################

## USEFUL FUNCTIONS WHEN CALLING make_basic_network ##

#####################################################################################################################

# !! These functions have to be called BEFORE calling make_basic_network !!

# map the index assigned to the bus in the basic network to the real one in the dataset
function basic_to_real_bus(data::Dict)

    bus_ordered = sort([bus for (i,bus) in data["bus"]], by=(x) -> x["index"])

    bus_id_map = Dict{Int,Int}()
    for (i,bus) in enumerate(bus_ordered)
        bus_id_map[i] = bus["index"]
    end

    return bus_id_map
end

function real_to_basic_bus(data::Dict)

    bus_ordered = sort([bus for (i,bus) in data["bus"]], by=(x) -> x["index"])

    bus_id_map = Dict{Int,Int}()
    for (i,bus) in enumerate(bus_ordered)
        bus_id_map[bus["index"]] = i 
    end

    return bus_id_map
end

# map the index assigned to the branch in the basic network to the real one in the dataset
function real_to_basic_branch(data::Dict)

    branch_ordered = sort([branch for (i,branch) in data["branch"]], by=(x) -> x["index"])

    branch_id_map = Dict{Int,Int}()
    for (i,bus) in enumerate(branch_ordered)
        branch_id_map[bus["index"]] = i 
    end

    return branch_id_map
end

function basic_to_real_branch(data::Dict)

    branch_ordered = sort([branch for (i,branch) in data["branch"]], by=(x) -> x["index"])

    branch_id_map = Dict{Int,Int}()
    for (i,bus) in enumerate(branch_ordered)
        branch_id_map[i] = bus["index"]
    end

    return branch_id_map
end


