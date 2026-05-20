using Pkg
install = 0
if install == 1
    Pkg.add("Parameters")
    Pkg.add("Plots")
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
    Pkg.add("QuantEcon")
    Pkg.add("DelimitedFiles")
end

using Parameters, Optim, Distributions, SharedArrays, Random, JLD2, Statistics, StatsBase, GLM, DataFrames, OrderedCollections, LinearAlgebra, FixedEffectModels, QuantEcon, DelimitedFiles
#Plots,  Distributed,
#addprocs(15)
include("mcmain_EC_model_annual_mod.jl")
#include("mcmain_EC_RFC_model_annual_mod.jl")
##############################################################
#####                     Optim                           ####
##############################################################
# x0_first3 = [1.0432552968181101, 2.5079538921320017, 0.7080712173179601, 0.2376451078412352] # not targeted
# opt_res_first3 = optimize(MSM_func_first3, x0_first3)
# minimizers_first3 = Optim.minimizer(opt_res_first3)
# minimizers_first3

# just for Optim, comment out below:


##############################################################
#####                     Optim Delta                     ####
##############################################################
rand_results = zeros(10, 6)
for i = 2:10
    runif = rand(Xoshiro(i), 5)
    prim, res = Initialize(3)
    println("Beginning of Iteration ", i, ":")
    # random_x0 = [runif[4] runif[5]*0.5 runif[1]*20+5 runif[2]*2 0.23 0.71 0.18]
    # random_x0 = [0.05 0.05 2.28655364976607 0.5038568304392312 0.15437472662956264 0.6878931557489516 0.1821507667406799]
    # random_x0 = [0.001 0.2*0.5038568304392312 2.28655364976607 0.5038568304392312] # 0.154067618167206 0.5385 0.05714972430845821]
    if i == 1
        random_x0 = [0.05903572698230247 0.02376585254741416 2.4327917664262735 0.5101173355824407] #0.154067618167206 0.5385 0.05714972430845821] 
    else
        random_x0 = [runif[1]*0.5 runif[2]*0.5 runif[3]*10+1 runif[4]] # 0.154067618167206 0.5385 0.05714972430845821]
    end
    opt_res_canon_random = optimize(MSM_delta_func_first3, random_x0)
    minimizers_canon_random = transpose(Optim.minimizer(opt_res_canon_random))
    println(minimizers_canon_random)
    open("export_capital.txt","a") do file
        println(file, minimizers_canon_random)
    end 
    rand_results[i,:] = hcat(Optim.minimum(opt_res_canon_random), minimizers_canon_random)
end
print(rand_results)

##############################################################
#####      Output Simulated Export Capital Panel          ####
##############################################################
# Initialize params and data
prim, res = Initialize(3)
panel = zeros(prim.n_sims*prim.n_firms*12, 6)
firm_export_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_labor_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_capital_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_sales_domestic = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_sales_all = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firms_export_sales = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
# Solve the model with the parameters
Solve_model(prim,res)
# Simulate the data
firm_export_choices, firm_labor_choices, firm_capital_choices, firm_sales_domestic, firm_sales_all, firms_export_sales = data_sim_delta_nsims(prim, res)
row = 1
for k = 1:prim.n_sims
    for j = 1:prim.n_firms
        for i = 1:12
            # Panel will be: Firm Year Export Capital Sales Export Sales
            panel[row,1] = parse(Float32, string(k,".",j,1))
            panel[row,2] = i
            panel[row,3] = firm_export_choices[100+i,j,k]
            panel[row,4] = firm_capital_choices[100+i,j,k]
            panel[row,5] = firm_sales_all[100+i,j,k]
            panel[row,6] = firms_export_sales[100+i,j,k]
            row += 1
        end
    end
end
writedlm("./data/Panel_Sim_delta_test.csv", panel, ",")
##############################################################
#####         Output Simulated Sunk Cost Panel            ####
##############################################################
prim, res = Initialize(1)
panel = zeros(prim.n_sims*prim.n_firms*12, 6)
firm_export_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_labor_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_capital_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_sales_domestic = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_sales_all = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
# Solve the model with the parameters
Solve_model(prim,res)
# Simulate the data
firm_export_choices, firm_labor_choices, firm_capital_choices, firm_sales_domestic, firm_sales_all, firms_export_sales = data_sim_delta_nsims(prim, res)
row = 1
for k = 1:prim.n_sims
    for j = 1:prim.n_firms
        for i = 1:12
            # Panel will be: Firm Year Export Capital Sales Export Sales
            panel[row,1] = parse(Float32, string(k,".",j,1))
            panel[row,2] = i
            panel[row,3] = firm_export_choices[100+i,j,k]
            panel[row,4] = firm_capital_choices[100+i,j,k]
            panel[row,5] = firm_sales_all[100+i,j,k]
            row += 1
        end
    end
end
writedlm("./data/Panel_Sim_noDelta.csv", panel, ",")
#######################################
# Export Concentration Export Capital:#
#######################################
# sum_dimension = 1 means across all years, = 2 means cross section
sum_dimension = 2
total_export_sales = sum(firms_export_sales[101:112,:,:], dims=sum_dimension)
if sum_dimension == 1
    export_concentration_p90 = zeros(prim.n_sims)
    for k = 1:prim.n_sims
        only_exporters_val = filter((x) -> x > 0, total_export_sales[1,:,k])
        total_exports_val = sum(total_export_sales[1,:,k])
        p90 = quantile!(only_exporters_val, 0.9)
        only_p90_exporters_val = filter((x) -> x >= p90, total_export_sales[1,:,k])
        total_exports_p90_val = sum(only_p90_exporters_val)
        export_concentration_p90[k] = total_exports_p90_val/total_exports_val
    end
    mean(export_concentration_p90)
else
    export_concentration_p90 = zeros(12, prim.n_sims)
    for k = 1:prim.n_sims
        for i = 1:12
            only_exporters_val = filter((x) -> x > 0, firms_export_sales[100+i,:,k])
            total_exports_val = total_export_sales[i,1,k]
            p90 = quantile!(only_exporters_val[:], 0.9)
            only_p90_exporters_val = filter((x) -> x >= p90, firms_export_sales[100+i,:,k])
            total_exports_p90_val = sum(only_p90_exporters_val[:])
            export_concentration_p90[i,k] = total_exports_p90_val/total_exports_val    
        end
    end
    mean(export_concentration_p90)
end
##################################
# Export Concentration Sunk Cost:#
##################################
# sum_dimension = 1 means across all years, = 2 means cross section
sum_dimension = 2
total_export_sales = sum(firms_export_sales[101:112,:,:], dims=sum_dimension)
if sum_dimension == 1
    export_concentration_p90 = zeros(prim.n_sims)
    for k = 1:prim.n_sims
        only_exporters_val = filter((x) -> x > 0, total_export_sales[1,:,k])
        total_exports_val = sum(total_export_sales[1,:,k])
        p90 = quantile!(only_exporters_val, 0.9)
        only_p90_exporters_val = filter((x) -> x >= p90, total_export_sales[1,:,k])
        total_exports_p90_val = sum(only_p90_exporters_val)
        export_concentration_p90[k] = total_exports_p90_val/total_exports_val
    end
    mean(export_concentration_p90)
else
    export_concentration_p90 = zeros(12, prim.n_sims)
    for k = 1:prim.n_sims
        for i = 1:12
            only_exporters_val = filter((x) -> x > 0, firms_export_sales[100+i,:,k])
            total_exports_val = total_export_sales[i,1,k]
            p90 = quantile!(only_exporters_val[:], 0.9)
            only_p90_exporters_val = filter((x) -> x >= p90, firms_export_sales[100+i,:,k])
            total_exports_p90_val = sum(only_p90_exporters_val[:])
            export_concentration_p90[i,k] = total_exports_p90_val/total_exports_val    
        end
    end
    mean(export_concentration_p90)
end