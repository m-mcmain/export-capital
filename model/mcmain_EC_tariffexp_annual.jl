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
# base = 1 # 1 for base model, 3 for delta
# filename = "base_annual"
# prim, res = Initialize(base) #initialize primitive and results structs for base model
# firms_export_capital_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_export_decisions_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_labor_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_capital_decisions_Base = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_sales_domestic_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_sales_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_export_sales_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# @time firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base, firms_export_sales_Base = tariff_experiment(prim, res, 10, 1, filename)
# Start_Base, Stop_Base, Foreign_Sales_Base, Sales_Base, Production_Base, output_Base = moment_calc_nsims(firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base)
# print(median(firms_sales_Base[prim.n_periods-7:prim.n_periods,:,:].-firms_sales_domestic_Base[prim.n_periods-7:prim.n_periods,:,:]))
# print("\\")

delta = 3
filename = "delta_annual"
prim, res = Initialize(delta) #initialize primitive and results structs for base model
firms_export_capital_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_decisions_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_labor_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_capital_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_domestic_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sales_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_export_sales_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_sunk_cost_spending_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
firms_profits_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
@time firms_export_decisions_Delta_noT, firms_labor_decisions_Delta_noT, firms_capital_decisions_Delta_noT, firms_sales_domestic_Delta_noT, firms_sales_Delta_noT, firms_export_sales_Delta_noT, firms_sunk_cost_spending_Delta_noT, firms_profits_Delta_noT, firms_export_capital_Delta_noT = tariff_experiment(prim, res, 0, 1, filename)
@time firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta, firms_export_sales_Delta, firms_sunk_cost_spending_Delta, firms_profits_Delta, firms_export_capital_Delta = tariff_experiment(prim, res, 10, 1, filename)
# Start_Delta, Stop_Delta, Foreign_Sales_Delta, Sales_Delta, Production_Delta, output_Delta = moment_calc_nsims(firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta)
# print(median(firms_sales_Base[prim.n_periods-7:prim.n_periods,:,:].-firms_sales_domestic_Base[prim.n_periods-7:prim.n_periods,:,:]))
# print("\\")

noT_exporting_shock = (firms_export_decisions_Delta_noT .- firms_export_decisions_Delta)[prim.n_periods-11,:,:] .> 0
noT_exporting_shock_any = (firms_export_decisions_Delta_noT .- firms_export_decisions_Delta)[prim.n_periods-11:prim.n_periods-8,:,:] .> 0
noT_exporting_shock_any_max = maximum(noT_exporting_shock_any,dims=1)[1,:,:]
sum(mean(noT_exporting_shock,dims=2))
mean(firms_export_decisions_Delta_noT[prim.n_periods-11,:,:] .> 0)
noT_exporting_preshock = firms_export_decisions_Delta_noT[prim.n_periods-12,:,:] .== 1 
T_exporting_preshock = firms_export_decisions_Delta[prim.n_periods-12,:,:] .== 1 
T_exporting_postshock = firms_export_decisions_Delta[prim.n_periods-7,:,:] .== 1
re_entrants_τ = noT_exporting_shock .* noT_exporting_preshock .* T_exporting_preshock .* T_exporting_postshock
sum(mean(re_entrants_τ,dims=2))
exiters_τ = noT_exporting_shock .* noT_exporting_preshock
sum(mean(exiters_τ,dims=2))
non_entrants_τ = BitMatrix(noT_exporting_shock .* (1 .- noT_exporting_preshock))
sum(mean(non_entrants_τ,dims=2))
out_after_τ = (firms_export_decisions_Delta_noT .- firms_export_decisions_Delta)[prim.n_periods-7,:,:] .> 0
median_firm_sales = median(firms_sales_Delta_noT[prim.n_periods-11,:,:])

# annual_sunk_cost_spending_diff = mean(sum(firms_sunk_cost_spending_Delta_noT[prim.n_periods-12:prim.n_periods,:,:], dims = 2) .- sum(firms_sunk_cost_spending_Delta[prim.n_periods-12:prim.n_periods,:,:], dims = 2), dims=3)
# average_annual_sunk_cost_spending_Delta = mean(sum(firms_sunk_cost_spending_Delta[prim.n_periods-12:prim.n_periods,:,:], dims = 2), dims=3)[1,1,1]
# periods_diff = range(-1,10,length=12)
# plot(periods_diff, annual_sunk_cost_spending_diff[1:12]./average_annual_sunk_cost_spending_Delta, xlabel="Year", ylabel="Scaled Difference in Sunk Cost Spending", linewidth=3, label="", dpi=300)
# savefig("./model/images/sunk_cost_compare_t_annual.png")
# periods_diff_lim = range(1,10,length=10)
# plot(periods_diff_lim, annual_sunk_cost_spending_diff[3:12]./average_annual_sunk_cost_spending_Delta, xlabel="Year", ylabel="Scaled Difference in Sunk Cost Spending", label="", ylim=(-0.55,0), linewidth=3, dpi=300)
# savefig("./model/images/sunk_cost_compare_t_annual_lim.png")

periods_diff = range(-1,10,length=12)
annual_profits_diff_all = sum(firms_profits_Delta_noT[prim.n_periods-12:prim.n_periods,noT_exporting_shock_any_max], dims = 2) .- sum(firms_profits_Delta[prim.n_periods-12:prim.n_periods,noT_exporting_shock_any_max], dims = 2)
plot(periods_diff, annual_profits_diff_all[1:12]./sum(firms_profits_Delta_noT[prim.n_periods-11:prim.n_periods,noT_exporting_shock_any_max], dims = 2), linewidth = 3, ylabel="Fraction of No Tax Profits", xlabel = "Year", label="", dpi=300)
savefig("./model/images/profits_compare_t_annual_all_4y.png")

annual_profits_diff_reentrants = sum(firms_profits_Delta_noT[prim.n_periods-12:prim.n_periods,re_entrants_τ], dims = 2) .- sum(firms_profits_Delta[prim.n_periods-12:prim.n_periods,re_entrants_τ], dims = 2)
plot(periods_diff, annual_profits_diff_reentrants[1:12]./sum(firms_profits_Delta_noT[prim.n_periods-11:prim.n_periods,re_entrants_τ], dims = 2), linewidth = 3, ylabel="Fraction of No Tax Profits", xlabel = "Year", label="", dpi=300)
savefig("./model/images/profits_compare_t_annual_reentrants_4y.png")

annual_sunk_cost_diff_reentrants = sum(firms_sunk_cost_spending_Delta[prim.n_periods-12:prim.n_periods,re_entrants_τ], dims = 2) .- sum(firms_sunk_cost_spending_Delta_noT[prim.n_periods-12:prim.n_periods,re_entrants_τ], dims = 2)
plot(periods_diff, annual_sunk_cost_diff_reentrants[1:12]./(1000*prim.w), linewidth = 3, ylabel="Labor Hours", xlabel = "Year", label="", dpi=300)
savefig("./model/images/sunk_cost_compare_t_annual_reentrants_4y.png")

annual_profits_diff_exiters = sum(firms_profits_Delta_noT[prim.n_periods-12:prim.n_periods,exiters_τ], dims = 2) .- sum(firms_profits_Delta[prim.n_periods-12:prim.n_periods,exiters_τ], dims = 2)
plot(periods_diff, annual_profits_diff_exiters[1:12]./sum(firms_profits_Delta_noT[prim.n_periods-11:prim.n_periods,exiters_τ], dims = 2), linewidth = 3, ylabel="Fraction of No Tax Profits", xlabel = "Year", label="", dpi=300)
savefig("./model/images/profits_compare_t_annual_exiters_4y.png")

annual_sunk_cost_diff_exiters = sum(firms_sunk_cost_spending_Delta[prim.n_periods-12:prim.n_periods,exiters_τ], dims = 2) .- sum(firms_sunk_cost_spending_Delta_noT[prim.n_periods-12:prim.n_periods,exiters_τ], dims = 2)
plot(periods_diff, annual_sunk_cost_diff_exiters[1:12]./(1000*prim.w), linewidth = 3, ylabel="Labor Hours", xlabel = "Year", label="", dpi=300)
savefig("./model/images/sunk_cost_compare_t_annual_exiters_4y.png")

annual_profits_diff_non_entrants = sum(firms_profits_Delta_noT[prim.n_periods-12:prim.n_periods,non_entrants_τ], dims = 2) .- sum(firms_profits_Delta[prim.n_periods-12:prim.n_periods,non_entrants_τ], dims = 2)
plot(periods_diff, annual_profits_diff_non_entrants[1:12]./sum(firms_profits_Delta_noT[prim.n_periods-11:prim.n_periods,non_entrants_τ], dims = 2), linewidth = 3, ylabel="Fraction of No Tax Profits", xlabel = "Year", label="", dpi=300)
savefig("./model/images/profits_compare_t_annual_non_entrants_4y.png")

annual_sunk_cost_subsidy = sum(res.prev_ex_grid[floor.(Int, firms_export_capital_Delta[prim.n_periods-7,out_after_τ])])/(1000*median_firm_sales)

periods_diff_lim = range(1,10,length=10)
plot(periods_diff_lim, annual_profits_diff[3:12], label="", dpi=300)
plot(periods_diff_lim, annual_profits_diff[3:12]./(average_annual_profits_Delta/2826), label="", dpi=300)

# periods = range(-1,5,length=7)
# plot(periods, [Start_Base[4:10] Start_Delta[4:10]], label=["Sunk Cost" "+ Export Capital"], dpi=300)
# xticks!(periods)
# savefig("./model/images/starter_compare_t_annual.png")
# plot(periods, [Stop_Base[1:12] Stop_Delta[1:12]], label=["Sunk Cost" "+ Export Capital"], dpi=300)
# xticks!(periods)
# savefig("./model/images/stopper_compare_t_annual.png")

# sales_Base_scaled = Sales_Base[5:11] ./ Sales_Base[5]
# sales_Delta_scaled = Sales_Delta[5:11] ./ Sales_Delta[5]
# sales_all = vcat(sales_Base_scaled, sales_Delta_scaled)
# plot(periods, [sales_Base_scaled sales_Delta_scaled], linewidth=3, label=["Sunk Cost" "+ Export Capital"], dpi=300)
# savefig("./model/images/sales_compare_t_annual.png")

# foreign_sales_Base_scaled = Foreign_Sales_Base[5:9] ./ Foreign_Sales_Base[5]
# foreign_sales_Delta_scaled = Foreign_Sales_Delta[5:9] ./ Foreign_Sales_Delta[5]
# foreign_sales_all = vcat(foreign_sales_Base_scaled, foreign_sales_Delta_scaled)
# groups = repeat(["δ=1", "Export Capital"], inner = 8)
# periods_full = repeat(-1:3)#, outer=2)
# plot(periods_full, [foreign_sales_Base_scaled foreign_sales_Delta_scaled], label=["Sunk Cost" "+ Export Capital"], dpi=300)
# savefig("./model/images/foreign_sales_compare_t_annual.png")

# firms_export_sales_Base_sumAvg = mean(sum(firms_export_sales_Base[prim.n_periods-11:prim.n_periods,:,:], dims = 2), dims=3)
# export_sales_Base_scaled = firms_export_sales_Base_sumAvg[5:12] ./ firms_export_sales_Base_sumAvg[5]

# firms_export_sales_Delta_sumAvg = mean(sum(firms_export_sales_Delta[prim.n_periods-11:prim.n_periods,:,:], dims = 2), dims=3)
# export_sales_Delta_scaled = firms_export_sales_Delta_sumAvg[5:12] ./ firms_export_sales_Delta_sumAvg[5]

# export_sales_all = vcat(export_sales_Base_scaled, export_sales_Delta_scaled)
# groups = repeat(["Sunk Cost", "+ Export Capital"], inner = 5)
# periods_full = repeat(-1:5)#, outer=2)
# periods_double = repeat(-1:5, outer=2)
# plot(periods_full, [export_sales_Base_scaled[1:7] export_sales_Delta_scaled[1:7]], linewidth=3, label=["Sunk Cost" "+ Export Capital"], ylabel="Scaled Export Revenues", xlabel="Period", dpi=300)
# savefig("./model/images/foreign_sales_compare_t_annual.png")

############################################################
#                    Q - Experiment                        #
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
# firms_export_sales_Base = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# @time firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base, firms_export_sales_Base = Q_experiment(prim, res, filename)
# Start_Base, Stop_Base, Foreign_Sales_Base, Sales_Base, Production_Base, output_Base = moment_calc_nsims(firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base)

# filename = "delta_annual"
# prim, res = Initialize(delta) #initialize primitive and results structs for base model
# firms_export_capital_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_export_decisions_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_labor_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_capital_decisions_Delta = ones(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_sales_domestic_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# firms_sales_Delta = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
# @time firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta = Q_experiment(prim, res, filename)
# Start_Delta, Stop_Delta, Foreign_Sales_Delta, Sales_Delta, Production_Delta, output_Delta = moment_calc_nsims(firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta)

# periods = range(-1,5,length=7)
# plot(periods, [Start_Base[4:10] Start_Delta[4:10]], linewidth=3, label=["Sunk Cost" "+ Export Capital"], dpi=300)
# xticks!(periods)
# savefig("./model/images/starter_compare_Q_annual.png")
# plot(periods, [Stop_Base[4:10] Stop_Delta[4:10]], linewidth=3, label=["Sunk Cost" "+ Export Capital"], dpi=300)
# xticks!(periods)
# savefig("./model/images/stopper_compare_Q_annual.png")

# sales_Base_scaled = Sales_Base[5:11] ./ Sales_Base[5]
# sales_Delta_scaled = Sales_Delta[5:11] ./ Sales_Delta[5]
# sales_all = vcat(sales_Base_scaled, sales_Delta_scaled)
# groups = repeat(["δ=1", "Export Capital"], inner = 8)
# periods_full = repeat(-1:5, outer=2)
# plot(periods, [sales_Base_scaled sales_Delta_scaled], linewidth=3, label=["Sunk Cost" "+ Export Capital"], dpi=300)
# savefig("./model/images/sales_compare_Q_annual.png")

# foreign_sales_Base_scaled = Foreign_Sales_Base[5:11] ./ Foreign_Sales_Base[5]
# foreign_sales_Delta_scaled = Foreign_Sales_Delta[5:11] ./ Foreign_Sales_Delta[5]
# foreign_sales_all = vcat(foreign_sales_Base_scaled, foreign_sales_Delta_scaled)
# groups = repeat(["δ=1", "Export Capital"], inner = 7)
# periods_full = repeat(-1:5)#, outer=2)
# periods_double = repeat(-1:5, outer=2)
# plot(periods_full, [foreign_sales_Base_scaled foreign_sales_Delta_scaled], linewidth=3, label=["Sunk Cost" "+ Export Capital"], ylabel="Scaled Export Revenues", dpi=300)
# savefig("./model/images/foreign_sales_compare_Q_annual.png")

###########################################################
#             Export Capital - Experiment                 #
###########################################################
base = 1 # 1 for base model, 3 for delta
delta = 3
filename = "base_annual"
prim, res = Initialize(base) #initialize primitive and results structs for base model
firms_export_capital_Base = ones(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_export_decisions_Base = zeros(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_labor_decisions_Base = ones(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_capital_decisions_Base = ones(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_sales_domestic_Base = zeros(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_sales_Base = zeros(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_export_sales_Base = zeros(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
@time firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base, firms_export_sales_Base = export_experience_experiment(prim, res, filename)
# Start_Base, Stop_Base, Foreign_Sales_Base, Sales_Base, Production_Base, output_Base = moment_calc_nsims(firms_export_decisions_Base, firms_labor_decisions_Base, firms_capital_decisions_Base, firms_sales_domestic_Base, firms_sales_Base)

val_annual_firms_export_decisions = zeros(112, prim.n_firms, prim.n_sims)
annual_firms_sales_foreign = zeros(112, prim.n_firms, prim.n_sims)
annual_total_sales_foreign = zeros(112)

for j = 1:112
    val_annual_firms_export_decisions[j,:,:] = val_annual_firms_export_decisions[j,:,:] .+ firms_export_decisions_Base[j+prim.n_periods-13,:,:]
    annual_firms_sales_foreign[j,:,:] = annual_firms_sales_foreign[j,:,:] .+ firms_sales_Base[j+prim.n_periods-13,:,:] .- firms_sales_domestic_Base[j+prim.n_periods-13,:,:]
end
 
annual_firms_export_decisions = val_annual_firms_export_decisions .> 0
mean_annual_firms_sales_foreign = mean(annual_firms_sales_foreign, dims=3)

for i = 1:112
    for j = 1:prim.n_firms
        if i == 1
            annual_total_sales_foreign[i] = annual_total_sales_foreign[i] + mean_annual_firms_sales_foreign[i,j]
        else
            annual_total_sales_foreign[i] = annual_total_sales_foreign[i] + mean_annual_firms_sales_foreign[i,j] 
        end
    end
end

foreign_sales_Base_scaled_long = annual_total_sales_foreign[1:112] ./ annual_total_sales_foreign[1]
periods_full = repeat(-1:13)
plot(periods_full, foreign_sales_Base_scaled_long[1:15], linewidth=3, label="", dpi=300)
savefig("./model/images/foreign_sales_return_EE_annual_base.png")


filename = "delta_annual"
prim, res = Initialize(delta) #initialize primitive and results structs for base model
firms_export_capital_Delta = ones(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_export_decisions_Delta = zeros(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_labor_decisions_Delta = ones(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_capital_decisions_Delta = ones(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_sales_domestic_Delta = zeros(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_sales_Delta = zeros(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
firms_export_sales_Delta = zeros(prim.n_periods_experiment+100, prim.n_firms, prim.n_sims)
@time firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta, firms_export_sales_Delta = export_experience_experiment(prim, res, filename)
# Start_Delta, Stop_Delta, Foreign_Sales_Delta, Sales_Delta, Production_Delta, output_Delta = moment_calc_nsims(firms_export_decisions_Delta, firms_labor_decisions_Delta, firms_capital_decisions_Delta, firms_sales_domestic_Delta, firms_sales_Delta)

val_annual_firms_export_decisions = zeros(112, prim.n_firms, prim.n_sims)
annual_firms_sales_foreign = zeros(112, prim.n_firms, prim.n_sims)
annual_total_sales_foreign = zeros(112)

for j = 1:112
    val_annual_firms_export_decisions[j,:,:] = val_annual_firms_export_decisions[j,:,:] .+ firms_export_decisions_Delta[j+prim.n_periods-13,:,:]
    annual_firms_sales_foreign[j,:,:] = annual_firms_sales_foreign[j,:,:] .+ firms_sales_Delta[j+prim.n_periods-13,:,:] .- firms_sales_domestic_Delta[j+prim.n_periods-13,:,:]
end
 
annual_firms_export_decisions = val_annual_firms_export_decisions .> 0
mean_annual_firms_sales_foreign = mean(annual_firms_sales_foreign, dims=3)

for i = 1:112
    for j = 1:prim.n_firms
        if i == 1
            annual_total_sales_foreign[i] = annual_total_sales_foreign[i] + mean_annual_firms_sales_foreign[i,j]
        else
            annual_total_sales_foreign[i] = annual_total_sales_foreign[i] + mean_annual_firms_sales_foreign[i,j] 
        end
    end
end

# periods = range(-1,5,length=7)
# plot(periods, [Start_Base[4:10] Start_Delta[4:10]], linewidth=3, label=["Sunk Cost" "+ Export Capital"], dpi=300)
# xticks!(periods)
# savefig("./model/images/starter_compare_EE_annual.png")
# plot(periods, [Stop_Base[4:10] Stop_Delta[4:10]], linewidth=3, label=["Sunk Cost" "+ Export Capital"], dpi=300)
# xticks!(periods)
# savefig("./model/images/stopper_compare_EE_annual.png")

# sales_Base_scaled = Sales_Base[5:11] ./ Sales_Base[5]
# sales_Delta_scaled = Sales_Delta[5:11] ./ Sales_Delta[5]
# sales_all = vcat(sales_Base_scaled, sales_Delta_scaled)
# groups = repeat(["δ=1", "Export Capital"], inner = 7)
# periods_full = repeat(-1:5, outer=2)
# plot(periods, [sales_Base_scaled sales_Delta_scaled], linewidth=3, label=["Sunk Cost" "+ Export Capital"], dpi=300)
# savefig("./model/images/sales_compare_EE_annual.png")

foreign_sales_Delta_scaled_long = annual_total_sales_foreign[1:85] ./ annual_total_sales_foreign[1]
periods_full = repeat(-1:83)
plot(periods_full, foreign_sales_Delta_scaled_long, linewidth=3, label="", dpi=300)
savefig("./model/images/foreign_sales_return_EE_annual_delta.png")

foreign_sales_Base_scaled = Foreign_Sales_Base[1:12] ./ Foreign_Sales_Base[1]
foreign_sales_Delta_scaled = Foreign_Sales_Delta[1:12] ./ Foreign_Sales_Delta[1]
foreign_sales_all = vcat(foreign_sales_Base_scaled, foreign_sales_Delta_scaled)
groups = repeat(["δ=1", "Export Capital"], inner = 7)
periods_full = repeat(-1:10)#, outer=2)
periods_double = repeat(-1:10, outer=2)
plot(periods_full, [foreign_sales_Base_scaled foreign_sales_Delta_scaled], linewidth=3, label=["Sunk Cost" "+ Export Capital"], dpi=300)
savefig("./model/images/foreign_sales_compare_EE_annual.png")

firms_export_sales_Base_sumAvg = mean(sum(firms_export_sales_Base[prim.n_periods-11:prim.n_periods,:,:], dims = 2), dims=3)
export_sales_Base_scaled = firms_export_sales_Base_sumAvg[5:11] ./ firms_export_sales_Base_sumAvg[5]

firms_export_sales_Delta_sumAvg = mean(sum(firms_export_sales_Delta[prim.n_periods-11:prim.n_periods,:,:], dims = 2), dims=3)
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
delta = 3

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

firms_export_sales_Delta_rand_sumAvg = mean(sum(firms_export_sales_Delta_rand[prim.n_periods-11:prim.n_periods,:,:], dims = 2), dims=3)
export_sales_rand_scaled = firms_export_sales_Delta_rand_sumAvg[5:11] ./ firms_export_sales_Delta_rand_sumAvg[5]

firms_export_sales_Delta_targeted_sumAvg = mean(sum(firms_export_sales_Delta_targeted[prim.n_periods-11:prim.n_periods,:,:], dims = 2), dims=3)
export_sales_targeted_scaled = firms_export_sales_Delta_targeted_sumAvg[5:11] ./ firms_export_sales_Delta_targeted_sumAvg[5]

export_sales_all = vcat(export_sales_rand_scaled, export_sales_targeted_scaled)
groups = repeat(["Random", "Targeted"], inner = 7)
periods_full = repeat(-1:5)
plot(periods_full, [export_sales_rand_scaled export_sales_targeted_scaled], linewidth=3, label=["Random" "Targeted"], ylabel="Scaled Export Revenues", xlabel="Period", dpi=300)
savefig("./model/images/foreign_sales_compare_subsidy_annual.png")

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