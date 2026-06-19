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
    Pkg.add("DelimitedFiles")
end

using Parameters, Optim, Distributions, SharedArrays, Random, JLD2, Statistics, StatsBase, GLM, DataFrames, OrderedCollections, LinearAlgebra, FixedEffectModels, QuantEcon, DelimitedFiles
#Plots,  Distributed,
#addprocs(15)
include("mcmain_EC_model_annual_mod.jl")
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
rand_results = zeros(10, 5)
for i = 13:20
    # Make sure these line up with model file
    model = 3
    model_file = "export_capital.txt"
    
    runif = rand(Xoshiro(i), 5)
    prim, res = Initialize(2)
    println("Beginning of Iteration ", i, ":")
    # random_x0 = [runif[4] runif[5]*0.5 runif[1]*20+5 runif[2]*2 0.23 0.71 0.18]
    random_x0 = [0.05 0.05 2.28655364976607 0.5038568304392312 0.15437472662956264 0.6878931557489516 0.1821507667406799]
    opt_res_canon_random = optimize(MSM_delta_func_first3, random_x0)
    minimizers_canon_random = transpose(Optim.minimizer(opt_res_canon_random))
    #println(minimizers_canon_random[1:3])
    println(Optim.minimum(opt_res_canon_random))
    rand_results[i-10,:] = vcat(Optim.minimum(opt_res_canon_random), minimizers_canon_random)
    open(model_file,"a") do file
        println(file, rand_results[i,:])
    end 
end
print(rand_results)

##############################################################
#####      Output Simulated Export Capital Panel          ####
##############################################################
# Initialize params and data
prim, res = Initialize(3)
panel = zeros(prim.n_sims*prim.n_firms*12, 7)
firm_export_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_labor_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_capital_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_sales_domestic = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_sales_all = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firms_export_sales = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
productivities = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
# Solve the model with the parameters
Solve_model(prim,res)
# Simulate the data
firm_export_choices, firm_labor_choices, firm_capital_choices, firm_sales_domestic, firm_sales_all, firms_export_sales, productivities = data_sim_delta_nsims_prod(prim, res)
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
            panel[row,7] = productivities[100+i,j,k]
            row += 1
        end
    end
end
writedlm("./data/Panel_Sim_delta.csv", panel, ",")
##############################################################
#####         Output Simulated Sunk Cost Panel            ####
##############################################################
prim, res = Initialize(1)
panel = zeros(prim.n_sims*prim.n_firms*12, 7)
firm_export_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_labor_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_capital_choices = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_sales_domestic = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
firm_sales_all = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
productivities = zeros(prim.n_periods, prim.n_firms, prim.n_sims)
# Solve the model with the parameters
Solve_model(prim,res)
# Simulate the data
firm_export_choices, firm_labor_choices, firm_capital_choices, firm_sales_domestic, firm_sales_all, firms_export_sales, productivities = data_sim_delta_nsims_prod(prim, res)
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
            panel[row,7] = productivities[100+i,j,k]
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