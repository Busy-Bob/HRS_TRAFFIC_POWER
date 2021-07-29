# Copyright (c) 2015: Changhyun Kwon.
# This file is part of https://github.com/chkwon/TrafficAssignment.jl/blob/master/src/load_network.jl 

# Traffic Assignment Data structure
mutable struct TA_Data
    network_name::String

    number_of_zones::Int
    number_of_nodes::Int
    first_thru_node::Int
    number_of_links::Int

    init_node::Array{Int,1}
    term_node::Array{Int,1}
    capacity::Array{Float64,1}
    link_length::Array{Float64,1}
    free_flow_time::Array{Float64,1}
    b::Array{Float64,1}
    power::Array{Float64,1}
    speed_limit::Array{Float64,1}
    toll::Array{Float64,1}
    link_type::Array{Int64,1}

    total_od_flow::Float64

    travel_demand::Array{Float64,2}
    od_pairs::Array{Tuple{Int64,Int64},1}

    toll_factor::Float64
    distance_factor::Float64

    best_objective::Float64
end

search_sc(s,c) = something(findfirst(isequal(c), s), 0)

function load_ta_network(network_name, network_data_file, trip_table_file; best_objective=-1.0, toll_factor=0.0, distance_factor=0.0)
    @assert ispath(network_data_file)
    @assert ispath(trip_table_file)

    ##################################################
    # Network Data
    ##################################################
    number_of_zones = 0
    number_of_links = 0
    number_of_nodes = 0
    first_thru_node = 0

    n = open(network_data_file, "r")

    while (line=readline(n)) != ""
        if occursin("<NUMBER OF ZONES>", line)
            number_of_zones = parse(Int, line[ search_sc(line, '>')+1 : end ] )
        elseif occursin("<NUMBER OF NODES>", line)
            number_of_nodes = parse(Int, line[ search_sc(line, '>')+1 : end ] )
        elseif occursin("<FIRST THRU NODE>", line)
            first_thru_node = parse(Int, line[ search_sc(line, '>')+1 : end ] )
        elseif occursin("<NUMBER OF LINKS>", line)
            number_of_links = parse(Int, line[ search_sc(line, '>')+1 : end ] )
        elseif occursin("<END OF METADATA>", line)
            break
        end
    end

    @assert number_of_links > 0

    init_node = Array{Int64}(undef, number_of_links)
    term_node = Array{Int64}(undef, number_of_links)
    capacity = zeros(number_of_links)
    link_length = zeros(number_of_links)
    free_flow_time = zeros(number_of_links)
    b = zeros(number_of_links)
    power = zeros(number_of_links)
    speed_limit = zeros(number_of_links)
    toll = zeros(number_of_links)
    link_type = Array{Int64}(undef, number_of_links)

    idx = 1
    while !eof(n)
      line = readline(n)
        if occursin("~", line) || line == ""
            continue
        end

        if occursin(";", line)
            line = strip(line, [' ', '\n', ';'])
            line = replace(line, ";" => "")

            numbers = split(line)
            init_node[idx] = parse(Int64, numbers[1])
            term_node[idx] = parse(Int64, numbers[2])
            capacity[idx] = parse(Float64, numbers[3])
            link_length[idx] = parse(Float64, numbers[4])
            free_flow_time[idx] = parse(Float64, numbers[5])
            b[idx] = parse(Float64, numbers[6])
            power[idx] = parse(Float64, numbers[7])
            speed_limit[idx] = parse(Float64, numbers[8])
            toll[idx] = parse(Float64, numbers[9])
            link_type[idx] = Int(round(parse(Float64, numbers[10])))

            idx = idx + 1
        end
    end

    ##################################################
    # Trip Table
    ##################################################

    number_of_zones_trip = 0
    total_od_flow = 0

    f = open(trip_table_file, "r")

    while (line=readline(f)) != ""
        if occursin("<NUMBER OF ZONES>", line)
            number_of_zones_trip = parse(Int, line[ search_sc(line, '>')+1 : end ] )
        elseif occursin("<TOTAL OD FLOW>", line)
            total_od_flow = parse(Float64, line[ search_sc(line, '>')+1 : end ] )
        elseif occursin("<END OF METADATA>", line)
            break
        end
    end

    @assert number_of_zones_trip == number_of_zones # Check if number_of_zone is same in both txt files
    @assert total_od_flow > 0

    travel_demand = zeros(number_of_zones, number_of_zones)
    od_pairs = Array{Tuple{Int64, Int64}}(undef, 0)

    origin = -1

    while !eof(f)
        line = readline(f)

        if line == ""
            origin = -1
            continue
        elseif occursin("Origin", line)
            origin = parse(Int, split(line)[2] )
        elseif occursin(";", line)
            pairs = split(line, ";")
            for i=1:size(pairs)[1]
                if occursin(":", pairs[i])
                    pair = split(pairs[i], ":")
                    destination = parse(Int64, strip(pair[1]) )
                    od_flow = parse(Float64, strip(pair[2]) )

                    # println("origin=$origin, destination=$destination, flow=$od_flow")

                    travel_demand[origin, destination] = od_flow
                    push!(od_pairs, (origin, destination))
                end
            end
        end
    end

    # Preparing data to return
    ta_data = TA_Data(
        network_name,
        number_of_zones,
        number_of_nodes,
        first_thru_node,
        number_of_links,
        init_node,
        term_node,
        capacity,
        link_length,
        free_flow_time,
        b,
        power,
        speed_limit,
        toll,
        link_type,
        total_od_flow,
        travel_demand,
        od_pairs,
        toll_factor,
        distance_factor,
        best_objective)

    return ta_data

end # end of load_network function



