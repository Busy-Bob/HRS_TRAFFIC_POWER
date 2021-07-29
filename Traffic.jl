module Traffic

using LightGraphs, SimpleWeightedGraphs
using Distributions
using JuMP
import Ipopt

export qʰALLTime, qᵍALLTime, solve_TA, Nₕ

##########################################
# SiouxFalls traffic info
##########################################
include("readTrafficInfo.jl")
const network_data_file = "./data/SiouxFalls_net.tntp";
const trip_table_file = "./data/SiouxFalls_trips.tntp";
td = load_ta_network("SiouxFalls", network_data_file, trip_table_file);
t₀ = td.free_flow_time;
cₐ = td.capacity;
n = td.number_of_nodes; 
b = td.number_of_links; 
roads = [start_node => end_node for (start_node, end_node) in zip(td.init_node, td.term_node)];
road2roadidx = Dict(roads[i]=>i  for i in 1:length(roads))

##########################################
# additional travel info
##########################################
BRP_int(xₐ, t₀, cₐ) = t₀ * (xₐ + 0.03 * xₐ^5  / cₐ^4) #  integral of BRP


# traffic graph
g = SimpleWeightedDiGraph(n);
for (start_node, end_node, link_length) in zip(td.init_node, td.term_node, td.link_length)
    add_edge!(g, start_node, end_node, link_length)
end
# Yen's algorithm to find k shortest paths
# - there is a bug in `rem_edge!` of SimpleWeightedGraphs.jl 
# - so convert SimpleWeightedGraphs to a SimpleDiGraph
# https://github.com/JuliaGraphs/LightGraphs.jl/issues/1505
# https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl/issues/66 
K = 10 # K shortest path in each OD
available_path_GV = Dict()
for start_node in 1:n
    for end_node in 1:n 
        if start_node ≠ end_node
            path_state = yen_k_shortest_paths(SimpleDiGraph(g), start_node, end_node, weights(g), K)
            temp_via_roads_idx = []
            for path in path_state.paths 
                via_roads = [path[i-1] => path[i] for i in 2:length(path)]
                via_roads_idx = [road2roadidx[road] for road in via_roads]
                push!(temp_via_roads_idx, via_roads_idx)
            end
            available_path_GV[start_node=>end_node] = temp_via_roads_idx
        end
    end
end

##########################################
# additional hydrogen refueling station info
##########################################
Eₕ = 10; # hydrogen required per vehicle: 10kg 
Nₕ = 8; # Number of hydrogen refueling stations

# ρⱼ = ones(Nₕ) * 10; # price of hydrogen kg/$
ω = 0.5 # price of time $/minute
#  hydrogen stations are on these roads

station_placement = [16=>8, 3=>4, 6=>5, 9=>10, 22=>21, 12=>11, 17=>19, 23=>14];
station2roadidx = [road2roadidx[road] for road in station_placement];

# find K paths that contain at least one HRS
available_path_FCEV = Dict()
available_path_via_stationidx_FCEV = Dict()
for start_node in 1:n
    for end_node in 1:n 
        if start_node ≠ end_node
            path_state = yen_k_shortest_paths(SimpleDiGraph(g), start_node, end_node, weights(g), 200)
            temp_via_roads_idx = []
            temp_via_stationidx = []
            # find K paths
            for path in path_state.paths 
                via_roads = [path[i-1] => path[i] for i in 2:length(path)]
                if (via_roads ∩ station_placement) ≠ []
                    via_roads_idx = [road2roadidx[road] for road in via_roads]
                    push!(temp_via_roads_idx, via_roads_idx)
                    push!(temp_via_stationidx, [findfirst(x->x==station, station_placement) 
                        for station in via_roads ∩ station_placement])
                    if length(temp_via_roads_idx) >= K
                        break
                    end
                end
            end
            available_path_FCEV[start_node=>end_node] = temp_via_roads_idx
            available_path_via_stationidx_FCEV[start_node=>end_node] = temp_via_stationidx
        end
    end
end

# from here, all dict keys are converted to the "index" in OD
# convenient to use in solve_TA
OD = [start_node=>end_node for start_node in 1:n for end_node in 1:n if start_node ≠ end_node]
# shape (OD_idx, path_idx, road)
pathʰ = [available_path_FCEV[OD_pairs] for OD_pairs in OD]
pathᵍ = [available_path_GV[OD_pairs] for OD_pairs in OD]
pathʰ_via_stationidx = [available_path_via_stationidx_FCEV[OD_pairs] for OD_pairs in OD]

# from road index to know all (OD, path) through it
road_idx2pathʰidx = [[] for _ in roads]
road_idx2pathᵍidx = [[] for _ in roads]
for (OD_idx, paths) in enumerate(pathʰ)
    for (path_idx, path) in enumerate(paths)
        for road_idx in path
            push!(road_idx2pathʰidx[road_idx], (OD=OD_idx, path=path_idx))
        end
    end
end
for (OD_idx, paths) in enumerate(pathᵍ)
    for (path_idx, path) in enumerate(paths)
        for road_idx in path
            push!(road_idx2pathᵍidx[road_idx], (OD=OD_idx, path=path_idx))
        end
    end
end

##########################################
# additional traffic demand info
##########################################
# assume the traffic demand is normal distribution
# assume 8:00 and 20:00 is the peak of the distribution
# its value: base_travel_demand * adjust_factor
# 2 a.m. to 2 p.m. : Normal(8, 3) * base_travel_demand * adjust_factor
# 2 p.m. to 2 a.m. : - Normal(20, 3) * base_travel_demand * adjust_factor
# the reverse traffic demand is the transpose of demand matrix
base_travel_demand = td.travel_demand
adjust_factor = 1
# travel demand for all FCEVs and GVs
travel_demand = [] # Vector{Matrix} (time, start_node, end_node)
T = 24
for t in 1:T
    if t < 2
        push!(travel_demand, pdf(Normal(20, 5), 24+t) / pdf(Normal(20, 5), 20) * base_travel_demand' * adjust_factor)
    elseif t >= 2 && t < 14
        push!(travel_demand, pdf(Normal(8, 5), t) / pdf(Normal(8, 5), 8)  * base_travel_demand * adjust_factor)
    else
        push!(travel_demand, pdf(Normal(20, 5), t) / pdf(Normal(20, 5), 20) * base_travel_demand' * adjust_factor)
    end
end
# FCEV retio, we think the FCEVs that don't require hydrogen as GVs
GV_ratio = 0.9998
FCEV_ration = 0.0002
# qʰALLTime shape: (T, OD)
qʰALLTime = [
    [ FCEV_ration * travel_demand[t][OD_pair.first, OD_pair.second] for OD_pair in OD] 
        for t in 1:T
]

qᵍALLTime = [ 
    [GV_ratio * travel_demand[t][OD_pair.first, OD_pair.second]  for OD_pair in OD]
        for t in 1:T
]


"""
- @description: solve a TA problem
- @param {ρⱼ} hydrogen price for each HRS
- @param {qᵍ} GV traffic demand for each OD
- @param {qʰ} FCEV traffic demand for each OD
- @return {hydrogen_demand}  FCEV refuling demand for each HRS
"""
function solve_TA(ρⱼ, qᵍ, qʰ)
    # pick the cheapest hydrogen station for each FCEV path
    pathʰ2stationidx = [ [ reduce((x, y) -> (ρⱼ[x] < ρⱼ[y]) ? x : y, stations_idx) 
        for stations_idx in paths] for paths in pathʰ_via_stationidx]

    # nonlinear optimization
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    register(model, :BRP_int, 3, BRP_int; autodiff = true)
    # JuMP.set_optimizer_attribute(model, "hessian_approximation", "limited-memory")


    ## decision variables
    @variable(model, xₐ[1:length(roads)] >= 0) # traffic flow for each road
    # traffic flow for each path
    fᵍ = [@variable(model, [1:length(pathᵍ[i])], lower_bound=0) for i in 1:length(OD)] # (OD, path) 
    fʰ = [@variable(model, [1:length(pathʰ[i])], lower_bound=0) for i in 1:length(OD)] # (OD, path)
    
    ## constraints
    # road flow is sum of path flows through the road
    @constraint(model, [a=1:length(roads)], 
        xₐ[a] == sum(fᵍ[idx.OD][idx.path] for idx in road_idx2pathᵍidx[a]) 
        + sum(fʰ[idx.OD][idx.path] for idx in road_idx2pathʰidx[a]))
    
    # OD demand is sum of their path flows
    @constraint(model, [i=1:length(OD)], qʰ[i] == sum(fʰ[i]))
    @constraint(model, [i=1:length(OD)], qᵍ[i] == sum(fᵍ[i]))

    ## objective
    @NLobjective(model, Min, 
        ω * sum(BRP_int(xₐ[i], t₀[i], cₐ[i]) for i in 1:length(roads)) # cost of time
        + Eₕ * sum( ρⱼ[station_idx] * fʰ[OD_idx][path_idx] # cost of hydrogen
            for (OD_idx, paths) in enumerate(pathʰ2stationidx)  
            for (path_idx, station_idx) in enumerate(paths) )
    )

    optimize!(model)
    println("【Traffic Result】:", raw_status(model))
    # if primal_status(model) == MOI.FEASIBLE_POINT
    #     println("Solution is a FEASIBLE POINT")
    # elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
    #     println("Solution is suboptimal due to a time limit, but a primal solution is available")
    #     println(ρⱼ)
    # else
    #     println("The model was not solved correctly.")
    #     println(ρⱼ)
    # end


    ## calculate the hydrogen demand for each HRS
    hydrogen_demand = zeros(Nₕ)
    for (OD_idx, paths) in enumerate(pathʰ2stationidx)  
        for (path_idx, station_idx) in enumerate(paths)
            hydrogen_demand[station_idx] += Eₕ * value(fʰ[OD_idx][path_idx])
        end
    end

    return hydrogen_demand
end

## test 
# solve_TA(ones(Nₕ), qᵍALLTime[1], qʰALLTime[1])

end


