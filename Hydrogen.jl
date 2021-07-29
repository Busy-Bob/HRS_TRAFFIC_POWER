module Hydrogen


using JuMP, Gurobi

export solve_HRS_op
#########################################
# hydrogen refueling station info
##########################################
T = 24; # time slot
Δt = 24 / T; 
Nₕ = 8; # the number of HRSs

p_max = 50 * ones(Nₕ) # Maximum hydrogen production power / MW
SOS₀ = 0.5 * ones(Nₕ) # initial State of storage
SOS_min = 0.05 * ones(Nₕ) 
SOS_max = 0.95 * ones(Nₕ) 

Wₛ = 50 * ones(Nₕ) # hydrogen tank / kg 
ηₕ = 0.95 # refueling efficiency (including compressor and refrigerator)
α = 1 / 60 * 1000 # conversion rate, kg H₂ / MWh
ρ₀ = 10 # average hydrogen price / $
ρ_max = 15
ρ_min = 5

β = 1e-5 # factor that balances hydrogen price gap and revenue
γ = 1e-9 # factor that smooths power 

# LMP and hydrogen demand obtained from other modules 
# shape (Nₕ, T)
# λⱼ = [ones(Nₕ, 8)  ones(Nₕ, 10)*1.6  ones(Nₕ, 6)*1.4]
# hydrogen_demand = [ones(Nₕ, 8)*30  ones(Nₕ, 10)*20  ones(Nₕ, 6)*30]

"""
- @description: solve an HRS operation problem
- @param {λⱼ} LMP price for each HRS dollar / MW
- @param {hydrogen_demand} hydrogen_demand for each HRS
- @return {ρⱼ}  hydrogen price for each HRS
- @return {pᵈ}  power for each HRS
"""
function solve_HRS_op(λⱼ, hydrogen_demand)

    model = Model(Gurobi.Optimizer)
    set_silent(model)

    @variable(model, ρⱼ[1:Nₕ, 1:T]); # hydrogen price
    @variable(model, pᵈ[1:Nₕ, 1:T]); # power from the grid
    @variable(model, Eₛ[1:Nₕ, 1:T]); # storage tank net increase

    @constraint(model, [j=1:Nₕ, t=1:T], Eₛ[j,t] == α * pᵈ[j,t] * Δt - hydrogen_demand[j,t] )
    @constraint(model, [j=1:Nₕ, τ=1:T], SOS₀[j] * Wₛ[j] + sum(Eₛ[j, 1:τ]) >= SOS_min[j] * Wₛ[j] )
    @constraint(model, [j=1:Nₕ, τ=1:T], SOS₀[j] * Wₛ[j] + sum(Eₛ[j, 1:τ]) <= SOS_max[j] * Wₛ[j] )
    @constraint(model, sum(Eₛ; dims=2) .== 0 )
    @constraint(model, 0 .<= pᵈ .<= p_max)
    @constraint(model, sum(ρⱼ; dims=1) .<= Nₕ * ρ₀)
    @constraint(model, ρ_min .<=  ρⱼ .<=  ρ_max )
    @constraint(model, sum(ρⱼ .* hydrogen_demand) <= sum(ρ₀ .* hydrogen_demand) )
    @objective(model, Max, sum(ρⱼ .* hydrogen_demand - λⱼ .* pᵈ .* Δt) 
        - β * sum((ρⱼ .- ρ₀) .* (ρⱼ .- ρ₀)) 
        - γ * sum(pᵈ .* pᵈ) 
    )

    optimize!(model)
    println("【Hydrogen Result】:", raw_status(model))
    if termination_status(model) == MOI.OPTIMAL
        println("Solution is optimal")
    elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
        println("Solution is suboptimal due to a time limit, but a primal solution is available")
        println(λⱼ)
        println(hydrogen_demand)
    else
        println("The model was not solved correctly.")
        println(λⱼ)
        println(hydrogen_demand)
    end

    display(value.(ρⱼ))
    display(value.(pᵈ))
    return (ρⱼ=value.(ρⱼ), pᵈ=value.(pᵈ))
end

## test 
# solve_HRS_op(λⱼ, hydrogen_demand)

end