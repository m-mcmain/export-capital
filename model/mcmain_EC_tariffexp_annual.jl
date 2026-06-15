using Pkg
install = 0
if install == 1
    Pkg.add("Parameters")
    Pkg.add("Plots")
    Pkg.add("StatsPlots")
    Pkg.add("Optim")
    Pkg.add("Distributions")
    Pkg.add("SharedArrays")
    Pkg.add("Distributed")
    Pkg.add("Random")
    Pkg.add("JLD2")
    Pkg.add("Statistics")
    Pkg.add("StatsBase")
    Pkg.add("GLM")
    Pkg.add("DataFrames")
    Pkg.add("OrderedCollections")
    Pkg.add("LinearAlgebra")
    Pkg.add("FixedEffectModels")
end
#cd(joinpath(pwd(), "model"))

using Parameters, Plots, StatsPlots, Optim, Distributions, SharedArrays, Distributed, Random, JLD2, Statistics, StatsBase, GLM, DataFrames, OrderedCollections, LinearAlgebra, FixedEffectModels, QuantEcon

#addprocs(15)
include("mcmain_EC_model_annual_mod.jl")
Random.seed!(17)

############################################################
#                 Tariff - Experiment                      #
# ############################################################
base = 1 # 1 for base model, 3 for delta
delta = 3
filename = "base_annual"
prim, res = Initialize(base) #initialize primitive and results structs for base model
firms_export_capital_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_decisions_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_labor_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_capital_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_domestic_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_sales_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
@time firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base, firms_export_sales_Base = tariff_experiment(prim, res, 10, 1, filename)
Start_Base, Stop_Base, Foreign_Sales_Base, Sales_Base, Production_Base, output_Base = moment_calc_nsims(firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base)
print(median(firms_sales_Base[105:112,:,:].-firms_sales_domestic_Base[105:112,:,:]))
print("\\")

filename = "delta_annual"
prim, res = Initialize(delta) #initialize primitive and results structs for base model
firms_export_capital_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_decisions_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_labor_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_capital_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_domestic_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_sales_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
@time firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta, firms_export_sales_Delta = tariff_experiment(prim, res, 10, 1, filename)
Start_Delta, Stop_Delta, Foreign_Sales_Delta, Sales_Delta, Production_Delta, output_Delta = moment_calc_nsims(firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta)
print(median(firms_sales_Base[105:112,:,:].-firms_sales_domestic_Base[105:112,:,:]))
print("\\")

periods = range(-1,5,length=7)
plot(periods, [Start_Base[4:10] Start_Delta[4:10].+0.01], label=["δ=1" "Export Capital"], ylabel="Starter Rate", xlabel="Period", dpi=300)
xticks!(periods)
savefig("./model/images/starter_compare_t_annual.png")
plot(periods, [Stop_Base[4:10] Stop_Delta[4:10]], label=["δ=1" "Export Capital"], ylabel="Starter Rate", xlabel="Period", dpi=300)
xticks!(periods)
savefig("./model/images/stopper_compare_t_annual.png")

sales_Base_scaled = Sales_Base[5:11] ./ Sales_Base[5]
sales_Delta_scaled = Sales_Delta[5:11] ./ Sales_Delta[5]
sales_all = vcat(sales_Base_scaled, sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 8)
periods_full = repeat(-1:6, outer=2)
plot(periods, [sales_Base_scaled sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./model/images/sales_compare_t_annual.png")

foreign_sales_Base_scaled = Foreign_Sales_Base[5:9] ./ Foreign_Sales_Base[5]
foreign_sales_Delta_scaled = Foreign_Sales_Delta[5:9] ./ Foreign_Sales_Delta[5]
foreign_sales_all = vcat(foreign_sales_Base_scaled, foreign_sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 8)
periods_full = repeat(-1:3)#, outer=2)
periods_double = repeat(-1:5, outer=2)
plot(periods_full, [foreign_sales_Base_scaled foreign_sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./model/images/foreign_sales_compare_t_annual.png")
# groupedbar(periods_double, foreign_sales_all, group = groups, label=["Export Capital" "δ=1"], dpi=300)
# savefig("./images/foreign_sales_compare_bar_t_annual.png")

# production_Base_scaled = Production_Base ./ maximum(Production_Base)
# production_Delta_scaled = Production_Delta ./ maximum(Production_Delta)
# production_all = vcat(production_Base_scaled, production_Delta_scaled)
# groups = repeat(["δ=1", "Export Capital"], inner = 8)
# periods_full = repeat(-1:6, outer=2)
# groupedbar(periods_full, production_all, group = groups, label=["δ=1" "Export Capital"], dpi=300)
# savefig("production_compare_annual.png")

firms_export_sales_Base_sumAvg = mean(sum(firms_export_sales_Base[101:112,:,:], dims = 2), dims=3)
export_sales_Base_scaled = firms_export_sales_Base_sumAvg[5:12] ./ firms_export_sales_Base_sumAvg[5]

firms_export_sales_Delta_sumAvg = mean(sum(firms_export_sales_Delta[101:112,:,:], dims = 2), dims=3)
export_sales_Delta_scaled = firms_export_sales_Delta_sumAvg[5:12] ./ firms_export_sales_Delta_sumAvg[5]

export_sales_all = vcat(export_sales_Base_scaled, export_sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 5)
periods_full = repeat(-1:6)#, outer=2)
periods_double = repeat(-1:5, outer=2)
plot(periods_full, [export_sales_Base_scaled export_sales_Delta_scaled], label=["δ=1" "Export Capital"], ylabel="Scaled Export Revenues", xlabel="Period", dpi=300)
savefig("./model/images/foreign_sales_compare_t_annual.png")

############################################################
#                    Q - Experiment                        #
############################################################

base = 1 # 1 for base model, 3 for delta
delta = 3
filename = "base_annual"
prim, res = Initialize(base) #initialize primitive and results structs for base model
firms_export_capital_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_decisions_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_labor_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_capital_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_domestic_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
@time firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base = Q_experiment(prim, res, filename)
Start_Base, Stop_Base, Foreign_Sales_Base, Sales_Base, Production_Base, output_Base = moment_calc_nsims(firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base)

filename = "delta_annual"
prim, res = Initialize(delta) #initialize primitive and results structs for base model
firms_export_capital_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_decisions_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_labor_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_capital_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_domestic_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
@time firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta = Q_experiment(prim, res, filename)
Start_Delta, Stop_Delta, Foreign_Sales_Delta, Sales_Delta, Production_Delta, output_Delta = moment_calc_nsims(firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta)

periods = range(-1,5,length=7)
plot(periods, [Start_Base[4:10] Start_Delta[4:10]], label=["δ=1" "Export Capital"], dpi=300)
xticks!(periods)
savefig("./model/images/starter_compare_Q_annual.png")
plot(periods, [Stop_Base[4:10] Stop_Delta[4:10]], label=["δ=1" "Export Capital"], dpi=300)
xticks!(periods)
savefig("./model/images/stopper_compare_Q_annual.png")

sales_Base_scaled = Sales_Base[5:11] ./ Sales_Base[5]
sales_Delta_scaled = Sales_Delta[5:11] ./ Sales_Delta[5]
sales_all = vcat(sales_Base_scaled, sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 8)
periods_full = repeat(-1:5, outer=2)
plot(periods, [sales_Base_scaled sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./model/images/sales_compare_Q_annual.png")
groupedbar(periods, [sales_Base_scaled sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./model/images/sales_compare_Q_annual_bar.png")

foreign_sales_Base_scaled = Foreign_Sales_Base[5:11] ./ Foreign_Sales_Base[5]
foreign_sales_Delta_scaled = Foreign_Sales_Delta[5:11] ./ Foreign_Sales_Delta[5]
foreign_sales_all = vcat(foreign_sales_Base_scaled, foreign_sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 7)
periods_full = repeat(-1:5)#, outer=2)
periods_double = repeat(-1:5, outer=2)
plot(periods_full, [foreign_sales_Base_scaled foreign_sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./model/images/foreign_sales_compare_Q_annual.png")
groupedbar(periods_double, foreign_sales_all, group = groups, label=["Export Capital" "δ=1"], dpi=300)
savefig("./model/images/foreign_sales_compare_bar_Q_annual.png")

############################################################
#              Export Capital - Experiment                 #
############################################################
base = 1 # 1 for base model, 3 for delta
delta = 3
filename = "base_annual"
prim, res = Initialize(base) #initialize primitive and results structs for base model
firms_export_capital_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_decisions_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_labor_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_capital_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_domestic_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_sales_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
@time firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base, firms_export_sales_Base = export_experience_experiment(prim, res, filename)
Start_Base, Stop_Base, Foreign_Sales_Base, Sales_Base, Production_Base, output_Base = moment_calc_nsims(firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base)

filename = "delta_annual"
prim, res = Initialize(delta) #initialize primitive and results structs for base model
firms_export_capital_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_decisions_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_labor_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_capital_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_domestic_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_sales_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
@time firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta, firms_export_sales_Delta = export_experience_experiment(prim, res, filename)
Start_Delta, Stop_Delta, Foreign_Sales_Delta, Sales_Delta, Production_Delta, output_Delta = moment_calc_nsims(firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta)

periods = range(-1,5,length=7)
plot(periods, [Start_Base[4:10] Start_Delta[4:10]], label=["δ=1" "Export Capital"], dpi=300)
xticks!(periods)
savefig("./model/images/starter_compare_EE_annual.png")
plot(periods, [Stop_Base[4:10] Stop_Delta[4:10]], label=["δ=1" "Export Capital"], dpi=300)
xticks!(periods)
savefig("./model/images/stopper_compare_EE_annual.png")

sales_Base_scaled = Sales_Base[5:11] ./ Sales_Base[5]
sales_Delta_scaled = Sales_Delta[5:11] ./ Sales_Delta[5]
sales_all = vcat(sales_Base_scaled, sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 7)
periods_full = repeat(-1:5, outer=2)
plot(periods, [sales_Base_scaled sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./model/images/sales_compare_EE_annual.png")
groupedbar(periods, [sales_Base_scaled sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./model/images/sales_compare_EE_annual_bar.png")

foreign_sales_Base_scaled = Foreign_Sales_Base[5:11] ./ Foreign_Sales_Base[5]
foreign_sales_Delta_scaled = Foreign_Sales_Delta[5:11] ./ Foreign_Sales_Delta[5]
foreign_sales_all = vcat(foreign_sales_Base_scaled, foreign_sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 7)
periods_full = repeat(-1:5)#, outer=2)
periods_double = repeat(-1:5, outer=2)
plot(periods_full, [foreign_sales_Base_scaled foreign_sales_Delta_scaled], label=["Sunk Cost" "+Export Capital"], dpi=300)
savefig("./model/images/foreign_sales_compare_EE_annual.png")
groupedbar(periods_double, foreign_sales_all, group = groups, label=["Export Capital" "δ=1"], dpi=300)
savefig("./model/images/foreign_sales_compare_bar_EE_annual.png")


firms_export_sales_Base_sumAvg = mean(sum(firms_export_sales_Base[101:112,:,:], dims = 2), dims=3)
export_sales_Base_scaled = firms_export_sales_Base_sumAvg[5:11] ./ firms_export_sales_Base_sumAvg[5]

firms_export_sales_Delta_sumAvg = mean(sum(firms_export_sales_Delta[101:112,:,:], dims = 2), dims=3)
export_sales_Delta_scaled = firms_export_sales_Delta_sumAvg[5:11] ./ firms_export_sales_Delta_sumAvg[5]

export_sales_all = vcat(export_sales_Base_scaled, export_sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 7)
periods_full = repeat(-1:5)#, outer=2)
periods_double = repeat(-1:5, outer=2)
plot(periods_full, [export_sales_Base_scaled export_sales_Delta_scaled], linewidth=3, label=["Sunk Cost" "+ Export Capital"], ylabel="Scaled Export Revenues", xlabel="Period", dpi=300)
savefig("./model/images/foreign_sales_compare_EE_annual.png")

############################################################
#              Variable Cost Uncertainty                   #
############################################################
# base = 1 # 1 for base model, 3 for delta
# delta = 3
# filename = "base_annual"
# prim, res = Initialize(base) #initialize primitive and results structs for base model
# firms_export_capital_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_export_decisions_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_labor_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_capital_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_sales_domestic_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_sales_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# @time firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base = tariff_uncertainty_experiment(prim, res, [0, 0.5], [[0.5, 0.5], [0.5, 0.5]], filename)
# Start_Base, Stop_Base, Foreign_Sales_Base, Sales_Base, Production_Base, output_Base = moment_calc_nsims(firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base)

# filename = "delta_annual"
# prim, res = Initialize(delta) #initialize primitive and results structs for base model
# firms_export_capital_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_export_decisions_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_labor_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_capital_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_sales_domestic_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_sales_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# @time firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta = tariff_uncertainty_experiment(prim, res, [0, 0.5], [[0.5, 0.5], [0.5, 0.5]], filename)
# Start_Delta, Stop_Delta, Foreign_Sales_Delta, Sales_Delta, Production_Delta, output_Delta = moment_calc_nsims(firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta)

# periods = range(-1,5,length=7)
# plot(periods, [Start_Base[4:10] Start_Delta[4:10]], label=["δ=1" "Export Capital"], dpi=300)
# xticks!(periods)
# savefig("./images/starter_compare_uncertain_annual.png")
# plot(periods, [Stop_Base[4:10] Stop_Delta[4:10]], label=["δ=1" "Export Capital"], dpi=300)
# xticks!(periods)
# savefig("./images/stopper_compare_uncertain_annual.png")

# sales_Base_scaled = Sales_Base[5:11] ./ Sales_Base[5]
# sales_Delta_scaled = Sales_Delta[5:11] ./ Sales_Delta[5]
# sales_all = vcat(sales_Base_scaled, sales_Delta_scaled)
# groups = repeat(["δ=1", "Export Capital"], inner = 7)
# periods_full = repeat(-1:5, outer=2)
# groupedbar(periods, [sales_Base_scaled sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
# savefig("./images/sales_compare_uncertain_annual.png")

# foreign_sales_Base_scaled = Foreign_Sales_Base[5:11] ./ Foreign_Sales_Base[5]
# foreign_sales_Delta_scaled = Foreign_Sales_Delta[5:11] ./ Foreign_Sales_Delta[5]
# foreign_sales_all = vcat(foreign_sales_Base_scaled, foreign_sales_Delta_scaled)
# groups = repeat(["δ=1", "Export Capital"], inner = 7)
# periods_full = repeat(-1:5)#, outer=2)
# periods_double = repeat(-1:5, outer=2)
# plot(periods_full, [foreign_sales_Base_scaled foreign_sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
# savefig("./images/foreign_sales_compare_uncertain_annual.png")
# groupedbar(periods_double, foreign_sales_all, group = groups, label=["Export Capital" "δ=1"], dpi=300)
# savefig("./images/foreign_sales_compare_bar_uncertain_annual.png")
############################################################
#              Export Subsidy Experiment                   #
############################################################
base = 1 # 1 for base model, 3 for delta
delta = 3
# prim, res = Initialize(base) #initialize primitive and results structs for base model
# firms_export_capital_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_export_decisions_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_labor_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_capital_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_sales_domestic_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_sales_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# @time firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base = export_subsidy_experiment(prim, res)

prim, res = Initialize(delta) #initialize primitive and results structs for base model

firms_export_capital_Delta_rand = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_decisions_Delta_rand = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_labor_decisions_Delta_rand = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_capital_decisions_Delta_rand = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_domestic_Delta_rand = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_Delta_rand = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_sales_Delta_rand = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)

firms_export_capital_Delta_targeted = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_decisions_Delta_targeted = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_labor_decisions_Delta_targeted = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_capital_decisions_Delta_targeted = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_domestic_Delta_targeted = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_Delta_targeted = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_sales_Delta_targeted = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)

productivity_and_exports = zeros(prim.n_firms, 3, prim.n_sims)

@time firms_export_decisions_Delta_rand, firms_labor_decisions_Delta_rand, firms_capital_decisions_Delta_rand, firms_sales_domestic_Delta_rand, firms_sales_Delta_rand, firms_export_sales_Delta_rand, firms_export_decisions_Delta_targeted, firms_labor_decisions_Delta_targeted, firms_capital_decisions_Delta_targeted, firms_sales_domestic_Delta_targeted, firms_sales_Delta_targeted, firms_export_sales_Delta_targeted, productivity_and_exports = export_subsidy_experiment(prim, res)

Start_Delta_rand, Stop_Delta_rand, Foreign_Sales_Delta_rand, Sales_Delta_rand, Production_Delta_rand, output_Delta_rand = moment_calc_nsims(firms_export_decisions_Delta_rand, firms_labor_decisions_Delta_rand, firms_capital_decisions_Delta_rand, firms_sales_domestic_Delta_rand, firms_sales_Delta_rand)
Start_Delta_targeted, Stop_Delta_targeted, Foreign_Sales_Delta_targeted, Sales_Delta_targeted, Production_Delta_targeted, output_Delta_targeted = moment_calc_nsims(firms_export_decisions_Delta_targeted, firms_labor_decisions_Delta_targeted, firms_capital_decisions_Delta_targeted, firms_sales_domestic_Delta_targeted, firms_sales_Delta_targeted)

periods = range(-1,5,length=7)
plot(periods, [Start_Delta_rand[4:10] Start_Delta_targeted[4:10]], label=["Random" "Targeted"], dpi=300)
xticks!(periods)
savefig("./model/images/starter_compare_subsidy_annual.png")
plot(periods, [Stop_Delta_rand[4:10] Stop_Delta_targeted[4:10]], label=["Random" "Targeted"], dpi=300)
xticks!(periods)
savefig("./model/images/stopper_compare_subsidy_annual.png")

periods = range(-1,5,length=7)
sales_rand_scaled = Sales_Delta_rand[5:11] ./ Sales_Delta_rand[5]
sales_targeted_scaled = Sales_Delta_targeted[5:11] ./ Sales_Delta_targeted[5]
sales_all = vcat(sales_rand_scaled, sales_targeted_scaled)
groups = repeat(["Random", "Targeted"], inner = 7)
periods_full = repeat(-1:5, outer=2)
plot(periods, [sales_rand_scaled sales_targeted_scaled], label=["Random" "Targeted"], dpi=300)
savefig("./model/images/sales_compare_subsidy_annual.png")
groupedbar(periods, [sales_rand_scaled sales_targeted_scaled], label=["Random" "Targeted"], dpi=300)
savefig("./model/images/sales_compare_subsidy_annual_bar.png")


firms_export_sales_Delta_rand_sumAvg = mean(sum(firms_export_sales_Delta_rand[101:112,:,:], dims = 2), dims=3)
export_sales_rand_scaled = firms_export_sales_Delta_rand_sumAvg[5:11] ./ firms_export_sales_Delta_rand_sumAvg[5]

firms_export_sales_Delta_targeted_sumAvg = mean(sum(firms_export_sales_Delta_targeted[101:112,:,:], dims = 2), dims=3)
export_sales_targeted_scaled = firms_export_sales_Delta_targeted_sumAvg[5:11] ./ firms_export_sales_Delta_targeted_sumAvg[5]

export_sales_all = vcat(export_sales_rand_scaled, export_sales_targeted_scaled)
groups = repeat(["Random", "Targeted"], inner = 7)
periods_full = repeat(-1:5)#, outer=2)
periods_double = repeat(-1:5, outer=2)
plot(periods_full, [export_sales_rand_scaled export_sales_targeted_scaled], linewidth=3, label=["Random" "Targeted"], ylabel="Scaled Export Revenues", xlabel="Period", dpi=300)
savefig("./model/images/foreign_sales_compare_subsidy_annual.png")
groupedbar(periods_double, export_sales_all, group = groups, label=["Random" "Targeted"], dpi=300)
#savefig("./model/images/foreign_sales_compare_bar_subsidy_annual.png")

##################
# ϵ-Distribution #
##################
exporters_rand_ind = productivity_and_exports[:,3,5] .== 1.0
exporters_prod_rand = productivity_and_exports[:,1,5][exporters_rand_ind]
non_exporters_rand_ind = productivity_and_exports[:,3,5] .== 0.0
non_exporters_prod_rand = productivity_and_exports[:,1,5][non_exporters_rand_ind]
b_range = range(0.3, 1.7, length = 15)
histogram(Any[exporters_prod_rand, non_exporters_prod_rand], fillcolor=[:blue :red], bins = b_range, fillalpha=0.4, dpi=300, label = ["" ""], xlabel = "Productivity", ylabel = "Count", title = "Random")
savefig("./model/images/subsidy_prod_dist_rand.png")

exporters_targeted_ind = productivity_and_exports[:,2,5] .== 1.0
exporters_prod_targeted = productivity_and_exports[:,1,5][exporters_targeted_ind]
non_exporters_targeted_ind = productivity_and_exports[:,2,5] .== 0.0
non_exporters_prod_targeted = productivity_and_exports[:,1,5][non_exporters_targeted_ind]
histogram(Any[exporters_prod_targeted, non_exporters_prod_targeted], fillcolor=[:blue :red], bins = b_range, fillalpha=0.4, dpi=300, label = ["Exporting" "Not Exporting"], xlabel = "Productivity", title="Targeted")
savefig("./model/images/subsidy_prod_dist_targeted.png")

histogram(Any[non_exporters_prod_rand, non_exporters_prod_targeted], ylims=(0,330), fillcolor=[:blue :red], bins = b_range, fillalpha=0.4, dpi=300, label = ["Random" "Targeted"], xlabel = "Productivity", ylabel = "Count", title="Firms Not Exporting")
savefig("./model/images/subsidy_prod_dist_non_exporters.png")
histogram(Any[exporters_prod_rand, exporters_prod_targeted], fillcolor=[:blue :red], bins = b_range, fillalpha=0.4, dpi=300, label = ["Random" "Targeted"], xlabel = "Productivity", ylabel = "Count", title="Firms Exporting")
savefig("./model/images/subsidy_prod_dist_exporters.png")

# Mean Additional Exporters:
mean(sum(productivity_and_exports[:,2,:] .- productivity_and_exports[:,3,:], dims=1))
# Mean Exporters:
mean(sum(productivity_and_exports[:,2,:], dims=1))