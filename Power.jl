module Power

using PowerSystems
using JuMP, Gurobi
using Random
using Distributions

export pᵈᵗAllTime, qᵈᵗAllTime, solve_LMP, bus, branch

Random.seed!(1234)
##########################################
# case22 info
##########################################
DATA_DIR = raw"./data/case22.m"
system_data = System(DATA_DIR)
# sorted bus and branch
bus = get_components(Bus, system_data) |> collect;
sort!(bus, by = x->x.number); # ascending order
branch = get_components(Branch, system_data) |> collect;
sort!(branch, by = branch->branch.arc.from.number);
# for each node, find the idx of the subsequent edge
π = [findall(x -> x.arc.from.number == i, branch)  for i in 1:length(bus)]
# voltage limit
Umin = [bus_i.voltage_limits.min^2 for bus_i in bus];
Umax = [bus_i.voltage_limits.max^2 for bus_i in bus];

##########################################
# case22 additional configuration
##########################################
# time slots
T = 24
# standard load ( mean = 1.0 )
standard_load = [0.80092226, 0.75151295, 0.72181773, 0.70716663, 0.7085614,
                0.72807473, 0.81076389, 0.90212761, 1.02712352, 1.12649447,
                1.17364172, 1.1766259 , 1.08926146, 1.08708576, 1.09357851,
                1.10167761, 1.13332739, 1.15798613, 1.15985351, 1.21728262,
                1.21873949, 1.16554919, 1.03838863, 0.90243688];
# MVA 
# In 2019, per-capita average energy consumption is 4.905 MWh in China
# 100,000 people city power load： 10^5 * 4.905 / 365 / 24 = 56 MW 
# power factor = 0.95, Q : P ≈ 0.3
bus_average_P_load = 10 * [ 0.000,  0.168,  0.168,  0.338,  0.146,  0.105,  0.088,  0.144, 
                 0.193,  0.144,  0.163,  0.163,  0.821,  0.347,  0.347,  0.803, 
                 0.496,  0.496,  0.438,  0.373,  0.373,  0.310];
bus_average_Q_load = 3 * [ 0.000,  0.209,  0.209,  0.373,  0.125,  0.142,  0.117,  0.186, 
                 0.259,  0.186,  0.195,  0.195,  0.717,  0.301,  0.301,  0.701, 
                 0.478,  0.478,  0.389,  0.360,  0.360,  0.294 ]; 
# traditional loads, (bus, T)
pᵈᵗAllTime =  (bus_average_P_load * standard_load') .* rand(Normal(1, 0.2), (length(bus), T))
qᵈᵗAllTime =  (bus_average_Q_load * standard_load') .* rand(Normal(1, 0.2), (length(bus), T))
# Generator limit / MVA 
pᵍmin = zeros(length(bus))
pᵍmax = ones(length(bus)) * 20
pᵍmax[1] = Inf # infinite input 
qᵍmin = ones(length(bus)) * -5
qᵍmax = ones(length(bus)) * 5
qᵍmax[1] = Inf # infinite input
# Generator cost coefficient 
a = 0.3 * ones(length(bus)) 
b = 150 * ones(length(bus)) 
# buying price of power market $ / MVA  
ρₘ = 140
# branch capacity limit / MVA 
Sˡ = ones(length(branch)) * 30;


"""
- @description: solve an LMP problem
- @param {pᵈᵗ} traditional active power demand for each bus
- @param {qᵈᵗ} traditional reactive power demand for each bus 
- @param {pᵈʰ} hydrogen electrolyser active power demand for each bus
- @param {qᵈʰ} hydrogen electrolyser reactive power demand for each bus
- @return {λ}  LMP price for each bus
"""
function solve_LMP(pᵈᵗ, qᵈᵗ, pᵈʰ,  qᵈʰ)
    # power load
    pᵈ = pᵈᵗ + pᵈʰ
    qᵈ = qᵈᵗ + qᵈʰ
    # optimization problem
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    JuMP.set_optimizer_attribute(model, "QCPDual", 1)
    # decision variables
    @variable(model, pᵍ[1:length(bus)]); #  generator active power for each bus
    @variable(model, qᵍ[1:length(bus)]); # generator reactive power for each bus
    @variable(model, U[1:length(bus)] >= 0); # squared voltage for each bus
    @variable(model, I[1:length(branch)] >= 0); # squared current for each branch
    @variable(model, Pˡ[1:length(branch)] >= 0); # active power of branch
    @variable(model, Qˡ[1:length(branch)] >= 0); # reactive power of branch
    # constraints
    @constraint(model, balance[k=1:length(branch)], Pˡ[k] + pᵍ[branch[k].arc.to.number] - branch[k].r * I[k] == sum(Pˡ[π[branch[k].arc.to.number]]) + pᵈ[branch[k].arc.to.number] )
    @constraint(model, [k=1:length(branch)], Qˡ[k] + qᵍ[branch[k].arc.to.number] - branch[k].x * I[k] == sum(Qˡ[π[branch[k].arc.to.number]]) + qᵈ[branch[k].arc.to.number] )
    @constraint(model, [k=1:length(branch)], U[branch[k].arc.to.number] == U[branch[k].arc.from.number] - 2 * (branch[k].r * Pˡ[k] + branch[k].x * Qˡ[k]) + (branch[k].r^2 + branch[k].x^2) * I[k] )
    @constraint(model, [k=1:length(branch)], [I[k] + U[branch[k].arc.from.number], 2Pˡ[k], 2Qˡ[k], I[k] - U[branch[k].arc.from.number]] ∈ SecondOrderCone())
    @constraint(model, pᵍmin .<= pᵍ .<= pᵍmax)
    @constraint(model, qᵍmin .<= qᵍ .<= qᵍmax)
    @constraint(model, Umin .<= U .<= Umax)
    @constraint(model, [k=1:length(branch)], [Sˡ[k], Pˡ[k], Qˡ[k]] ∈ SecondOrderCone())
    @constraint(model, [k=1:length(branch)], Pˡ[k] - branch[k].r * I[k] >= 0 )
    @constraint(model, [k=1:length(branch)], Qˡ[k] - branch[k].x * I[k] >= 0 )
    # objective
    @objective(model, Min, sum(a .* pᵍ .^ 2 + b .* pᵍ) + ρₘ * sum(Pˡ[π[1]]) )

    optimize!(model)
    println("【Power Result】:", raw_status(model))
    # if termination_status(model) == MOI.OPTIMAL
    #     println("Solution is optimal")
    # elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
    #     println("Solution is suboptimal due to a time limit, but a primal solution is available")
    #     println(pᵈʰ)
    # else
    #     println("The model was not solved correctly.")
    #     println(pᵈʰ)
    # end

    return dual.(balance)
end


end
