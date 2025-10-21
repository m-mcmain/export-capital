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
using Parameters, Plots, StatsPlots, Optim, Distributions, SharedArrays, Distributed, Random, JLD2, Statistics, StatsBase, GLM, DataFrames, OrderedCollections, LinearAlgebra, FixedEffectModels, QuantEcon

#addprocs(15)
include("mcmain_EC_model_annual_mod.jl")
Random.seed!(17)
############################################################
#                 Tariff - Experiment                      #
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
@time firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base = tariff_experiment(prim, res, 10, 1, filename)
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
@time firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta = tariff_experiment(prim, res, 10, 0, filename)
Start_Delta, Stop_Delta, Foreign_Sales_Delta, Sales_Delta, Production_Delta, output_Delta = moment_calc_nsims(firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta)
print(median(firms_sales_Base[105:112,:,:].-firms_sales_domestic_Base[105:112,:,:]))
print("\\")

periods = range(-1,6,length=8)
plot(periods, [Start_Base[4:11] Start_Delta[4:11]], label=["δ=1" "Export Capital"], dpi=300)
xticks!(periods)
savefig("./images/starter_compare_t_annual.png")
plot(periods, [Stop_Base[4:11] Stop_Delta[4:11]], label=["δ=1" "Export Capital"], dpi=300)
xticks!(periods)
savefig("./images/stopper_compare_t_annual.png")

sales_Base_scaled = Sales_Base[5:12] ./ Sales_Base[5]
sales_Delta_scaled = Sales_Delta[5:12] ./ Sales_Delta[5]
sales_all = vcat(sales_Base_scaled, sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 8)
periods_full = repeat(-1:6, outer=2)
plot(periods, [sales_Base_scaled sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./images/sales_compare_t_annual.png")

foreign_sales_Base_scaled = Foreign_Sales_Base[5:12] ./ Foreign_Sales_Base[5]
foreign_sales_Delta_scaled = Foreign_Sales_Delta[5:12] ./ Foreign_Sales_Delta[5]
foreign_sales_all = vcat(foreign_sales_Base_scaled, foreign_sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 8)
periods_full = repeat(-1:6)#, outer=2)
periods_double = repeat(-1:6, outer=2)
plot(periods_full, [foreign_sales_Base_scaled foreign_sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./images/foreign_sales_compare_t_annual.png")
groupedbar(periods_double, foreign_sales_all, group = groups, label=["Export Capital" "δ=1"], dpi=300)
savefig("./images/foreign_sales_compare_bar_t_annual.png")
#=
production_Base_scaled = Production_Base ./ maximum(Production_Base)
production_Delta_scaled = Production_Delta ./ maximum(Production_Delta)
production_all = vcat(production_Base_scaled, production_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 8)
periods_full = repeat(-1:6, outer=2)
groupedbar(periods_full, production_all, group = groups, label=["δ=1" "Export Capital"], dpi=300)
savefig("production_compare_annual.png")
=#
############################################################
#                    Q - Experiment                        #
############################################################
#=
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

periods = range(-1,6,length=8)
plot(periods, [Start_Base[4:11] Start_Delta[4:11]], label=["δ=1" "Export Capital"], dpi=300)
xticks!(periods)
savefig("./images/starter_compare_Q_annual.png")
plot(periods, [Stop_Base[4:11] Stop_Delta[4:11]], label=["δ=1" "Export Capital"], dpi=300)
xticks!(periods)
savefig("./images/stopper_compare_Q_annual.png")

sales_Base_scaled = Sales_Base[5:12] ./ Sales_Base[5]
sales_Delta_scaled = Sales_Delta[5:12] ./ Sales_Delta[5]
sales_all = vcat(sales_Base_scaled, sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 8)
periods_full = repeat(-1:6, outer=2)
groupedbar(periods, [sales_Base_scaled sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./images/sales_compare_Q_annual.png")

foreign_sales_Base_scaled = Foreign_Sales_Base[5:12] ./ Foreign_Sales_Base[5]
foreign_sales_Delta_scaled = Foreign_Sales_Delta[5:12] ./ Foreign_Sales_Delta[5]
foreign_sales_all = vcat(foreign_sales_Base_scaled, foreign_sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 8)
periods_full = repeat(-1:6)#, outer=2)
periods_double = repeat(-1:6, outer=2)
plot(periods_full, [foreign_sales_Base_scaled foreign_sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./images/foreign_sales_compare_Q_annual.png")
groupedbar(periods_double, foreign_sales_all, group = groups, label=["Export Capital" "δ=1"], dpi=300)
savefig("./images/foreign_sales_compare_bar_Q_annual.png")
=#
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
@time firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base = export_experience_experiment(prim, res, filename)
Start_Base, Stop_Base, Foreign_Sales_Base, Sales_Base, Production_Base, output_Base = moment_calc_nsims(firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base)

filename = "delta_annual"
prim, res = Initialize(delta) #initialize primitive and results structs for base model
firms_export_capital_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_decisions_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_labor_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_capital_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_domestic_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
@time firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta = export_experience_experiment(prim, res, filename)
Start_Delta, Stop_Delta, Foreign_Sales_Delta, Sales_Delta, Production_Delta, output_Delta = moment_calc_nsims(firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta)

periods = range(-1,6,length=8)
plot(periods, [Start_Base[4:11] Start_Delta[4:11]], label=["δ=1" "Export Capital"], dpi=300)
xticks!(periods)
savefig("./images/starter_compare_EE_annual.png")
plot(periods, [Stop_Base[4:11] Stop_Delta[4:11]], label=["δ=1" "Export Capital"], dpi=300)
xticks!(periods)
savefig("./images/stopper_compare_EE_annual.png")

sales_Base_scaled = Sales_Base[5:12] ./ Sales_Base[5]
sales_Delta_scaled = Sales_Delta[5:12] ./ Sales_Delta[5]
sales_all = vcat(sales_Base_scaled, sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 8)
periods_full = repeat(-1:6, outer=2)
groupedbar(periods, [sales_Base_scaled sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./images/sales_compare_EE_annual.png")

foreign_sales_Base_scaled = Foreign_Sales_Base[5:12] ./ Foreign_Sales_Base[5]
foreign_sales_Delta_scaled = Foreign_Sales_Delta[5:12] ./ Foreign_Sales_Delta[5]
foreign_sales_all = vcat(foreign_sales_Base_scaled, foreign_sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 8)
periods_full = repeat(-1:6)#, outer=2)
periods_double = repeat(-1:6, outer=2)
plot(periods_full, [foreign_sales_Base_scaled foreign_sales_Delta_scaled], label=["δ=1" "Export Capital"], dpi=300)
savefig("./images/foreign_sales_compare_EE_annual.png")
groupedbar(periods_double, foreign_sales_all, group = groups, label=["Export Capital" "δ=1"], dpi=300)
savefig("./images/foreign_sales_compare_bar_EE_annual.png")