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
    Pkg.add("MacroEconometricModels")
    Pkg.add("StatFiles")
    Pkg.add("MixedModels")
end

using Plots, Parameters, Optim, Distributions, SharedArrays, Random, JLD2, Statistics, StatsBase, GLM, DataFrames, OrderedCollections, LinearAlgebra, FixedEffectModels, QuantEcon, DelimitedFiles, StatFiles, MixedModels
#, MacroEconometricModels, Distributed,
#addprocs(15)
include("mcmain_EC_model_annual_mod.jl")


data_df = DataFrame(load("/home/m/mcmain/export-capital/data/dta/Colombia_EAM_panel_regReady_backup_julia.dta"))

pd_data = xtset(data_df, :nordest, :year)
m_probit_data = estimate_xtprobit(pd_data, :exported_est, [:exported_est_prev, :exported_est_prev2, :exported_est_prev3, 
                                                           :exported_est_prev4, :exported_est_prev5, :exported_est_prev6,
                                                           :exported_est_prev7, :log_sales_est_lag, :log_average_wage_lag, :log_total_capital_lag,
                                                           :exported_est_first_max, :mean_log_sales_est, :log_sales_est, :mean_log_average_wage_lag,
                                                           :mean_log_total_capital_lag, :mean_log_sales_est_lag, :y_2016, :y_2017, :y_2018, :y_2019,
                                                           :ind_10, :ind_11, :ind_13, :ind_14, :ind_15, :ind_16, :ind_17, :ind_18, :ind_19, :ind_20,
                                                           :ind_21, :ind_22, :ind_23, :ind_24, :ind_25, :ind_27, :ind_28, :ind_29, :ind_31, :ind_32]; model=:re, maxiter=1000)
m_ame_data = marginal_effects(m_probit_data) 
MacroEconometricModels.report(m_probit_data)
MacroEconometricModels.report(m_ame_data)

# Try with GLM
form = @formula(exported_est ~ 1 + exported_est_prev + exported_est_prev2 + exported_est_prev3 + 
                                exported_est_prev4 + exported_est_prev5 + exported_est_prev6 +
                                exported_est_prev7 + log_sales_est_lag + log_average_wage_lag + log_total_capital_lag +
                                exported_est_first_max + mean_log_sales_est + log_sales_est + mean_log_average_wage_lag +
                                mean_log_total_capital_lag + mean_log_sales_est_lag + y_2016 + y_2017 + y_2018 + y_2019 +
                                ind_10 + ind_11 + ind_13 + ind_14 + ind_15 + ind_16 + ind_17 + ind_18 + ind_19 + ind_20 +
                                ind_21 + ind_22 + ind_23 + ind_24 + ind_25 + ind_27 + ind_28 + ind_29 + ind_31 + ind_32 + (1 | nordest))

# 2. Fit the model using a Binomial distribution and a Probit link function
model = fit(GeneralizedLinearMixedModel, form, data_df, Binomial(), ProbitLink())

# 3. Display the summary results
println(model)

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
for i = 2:5
    # Make sure these line up with model file
    model = 3
    model_file = "export_capital.txt"
    
    runif = rand(Xoshiro(i), 5)
    prim, res = Initialize(model)
    println("Beginning of Iteration ", i, ":")
    # random_x0 = [runif[4] runif[5]*0.5 runif[1]*20+5 runif[2]*2 0.23 0.71 0.18]
    random_x0 = [0.03897024180506729 0.007706900420490436 2.8967193278506573 0.5160469782313699]
    # random_x0 = [runif[1]*0.2 runif[2]*0.1 runif[3]*5 runif[4]]
    opt_res_canon_random = optimize(MSM_delta_func_first3, random_x0)
    minimizers_canon_random = transpose(Optim.minimizer(opt_res_canon_random))
    #println(minimizers_canon_random[1:3])
    println(Optim.minimum(opt_res_canon_random))
    rand_results[i,:] = vcat(Optim.minimum(opt_res_canon_random), minimizers_canon_random)
    open(model_file,"a") do file
        println(file, rand_results[i,:])
    end 
    rand_results[i] = hcat(Optim.minimum(opt_res_canon_random), minimizers_canon_random)
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
firm_export_choices, firm_labor_choices, firm_capital_choices, firm_sales_domestic, firm_sales_all, firms_export_sales, productivities = data_sim_delta_nsims(prim, res)

###
# Create %Δ Plots for Delta
###
lchange_sales_Delta = zeros(prim.n_periods_experiment-1, prim.n_firms, prim.n_sims)
lchange_sales_Delta = (firm_sales_all[2:prim.n_periods_experiment,:,:].-firm_sales_all[1:prim.n_periods_experiment-1,:,:])./firm_sales_all[1:prim.n_periods_experiment-1,:,:]
avg_lchange_sales_Delta = mean(lchange_sales_Delta, dims = [2,3])
plot(avg_lchange_sales_Delta[2:prim.n_periods_experiment-1])

total_export_sales_Delta = sum(firms_export_sales, dims = 2)
plot(mean(total_export_sales_Delta,dims=3)[:,1,1])

total_sales_Delta = sum(firm_sales_all, dims = 2)
plot(mean(total_sales_Delta,dims=3)[:,1,1])

lchange_export_sales_Delta = zeros(prim.n_periods_experiment-1, prim.n_firms, prim.n_sims)
lchange_export_sales_Delta = (firms_export_sales[2:prim.n_periods_experiment,:,:].-firms_export_sales[1:prim.n_periods_experiment-1,:,:])./firms_export_sales[1:prim.n_periods_experiment-1,:,:]
avg_lchange_export_sales_Delta = mean(lchange_export_sales_Delta, dims = [2,3])
plot(avg_lchange_export_sales_Delta[2:prim.n_periods_experiment-1])

## Starter/Stopper Rate Over Time
firms_export_decisions_no_zero = hcat(firm_export_choices, ones(prim.n_periods_experiment, 1, prim.n_sims), zeros(prim.n_periods_experiment, 1, prim.n_sims)) 
firms_export_change = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
for i = 2:prim.n_periods_experiment
    for j = 1:prim.n_firms
        firms_export_change[i,j,:] = firm_export_choices[i,j,:] .- firm_export_choices[i-1,j,:]
    end
end
firms_export_change_no_zero = hcat(firms_export_change, ones(prim.n_periods_experiment, 1, prim.n_sims), -1 .* ones(prim.n_periods_experiment, 1, prim.n_sims)) 
starter_rate = zeros(prim.n_periods_experiment-1,prim.n_sims)
stopper_rate = zeros(prim.n_periods_experiment-1,prim.n_sims)

exporter_i = countmap(firms_export_change_no_zero[1,:,:])
exporter_i_minus1 = countmap(firms_export_change_no_zero[1,:,:])
for i = 2:prim.n_periods_experiment-1
    for k = 1:prim.n_sims
        exporter_i_minus1 = countmap(firms_export_decisions_no_zero[i-1,:,k])
        exporter_i = countmap(firms_export_change_no_zero[i,:,k])
        
        starter_rate[i-1,k] = exporter_i[1]/exporter_i_minus1[0]
        stopper_rate[i-1,k] = exporter_i[-1]/exporter_i_minus1[1]
    end
end
avg_starter_rate_Delta = mean(starter_rate, dims=2)
plot(avg_starter_rate_Delta)
avg_stopper_rate_Delta = mean(stopper_rate,dims=2)
plot(avg_stopper_rate_Delta)

## Coef of Variation Over Time
coef_variation = zeros(prim.n_periods_experiment-20, prim.n_sims)
for j = 1:prim.n_sims
    for i = 1:prim.n_periods_experiment-20
        coef_variation[i,j] = std(log.(firm_sales_domestic[i:i+20,:,j]))/mean(log.(firm_sales_domestic[i:i+20,:,j]))
    end
end
coef_variation_delta = mean(coef_variation,dims=2)
plot(coef_variation_delta)

## β-Coef Over Time
β_sales_delta = zeros(prim.n_periods_experiment-20, 1)
for i = 2:prim.n_periods_experiment-20
    for s = 1:1
        sales_today = firm_sales_domestic[i+1,:,s]
        sales_yesterday = firm_sales_domestic[i,:,s]
        year = []
        firm_num = []
        for i = i+2:i+19
            sales_today = [sales_today; firm_sales_domestic[i,:,s]]
            sales_yesterday = [sales_yesterday; firm_sales_domestic[i-1,:,s]]
        end

        for i = 2:20
            for j = 1:prim.n_firms
                year = [year; string("Y",i)]
                firm_num = [firm_num; string("F",j)]
            end
        end

        sales_today_vec = copyto!(Vector{Float64}(undef,length(sales_today)),sales_today)
        sales_yesterday_vec = copyto!(Vector{Float64}(undef,length(sales_yesterday)),sales_yesterday)

        reg_DF = DataFrame(s_today = log.(sales_today_vec), s_yesterday = log.(sales_yesterday_vec), y = year, f = firm_num)
        reg_results = reg(reg_DF, @formula(s_today ~ s_yesterday + fe(y) + fe(f)))

        β_moment = coef(reg_results)
        β_sales_delta[i,s] = β_moment[1]
    end
end
plot(β_sales_delta)

## Export-Sales Ratio Over Time
export_sales_ratio_delta = zeros(prim.n_periods_experiment-20)
for i = 1:prim.n_periods_experiment-20
    export_sales_ratio_delta[i] = sum(firms_export_sales[i:i+19,:,:]./firm_sales_all[i:i+19,:,:])/sum(firm_export_choices[i:i+19,:,:])
end
plot(export_sales_ratio_delta)

## Average Years Out and Frac Immedate Re-Entry Over Time
avg_years_out_delta = zeros(prim.n_periods_experiment-12)
frac_immediate_reentry_delta = zeros(prim.n_periods_experiment-12)
reentries = zeros(prim.n_periods_experiment-12)
years_out_before_reentry = zeros(prim.n_periods_experiment-12)
exits_all_noEndYears = zeros(prim.n_periods_experiment-12)
immediate_reentries = zeros(prim.n_periods_experiment-12)

for i=1:prim.n_periods_experiment-12
    
    exits_all = Int.(firms_export_change[i:i+11,:,:] .== -1)
    exits_cumulative = Int.(firms_export_change[i:i+11,:,:] .== -1)
    final_exit_no_reentry = ones(Int, 12, prim.n_firms, prim.n_sims)
    entries_all = Int.(firms_export_change[i:i+11,:,:] .== 1)

    for j in 2:12
        exits_cumulative[j,:,:] .+= exits_cumulative[j-1,:,:]
    end

    # Set entries where the firm is exporting equal to zero, keeping strings of numbers with length equal to years out of each exit
    exits_cumulative = exits_cumulative .* (ones(Int, 12, prim.n_firms, prim.n_sims) .- firm_export_choices[i:i+11,:,:])
    # Test if that exit is the final exit and has no re-entry tied to it
    for j in 1:11
        final_exit_no_reentry[j,:,:] = Int.((exits_cumulative[j,:,:] .== exits_cumulative[12,:,:]))
    end
    # Only keep exits that result in re-entries
    exits_cumulative = exits_cumulative .* (1 .- final_exit_no_reentry)
    # Loop through the number of exits and sum up the years out that result in re-entry
    max_exits = maximum(exits_cumulative)
    exits_count_map = countmap(exits_cumulative)
    years_out_before_reentry[i] = 0
    for j in 1:Int(max_exits)
        # Test if that number of exits actually occurs ever
        try
            exits_count_map[j]
        catch
            # if not set it to 0
            exits_count_map[j] = 0
        end
        years_out_before_reentry[i] += exits_count_map[j]
    end
    # Get maximum exits for each firm-simulation, which is number of re-entries
    reentries[i] = sum(maximum(exits_cumulative, dims = 1))
    # If there are no re-entries make them very small so it will be very far off the moment
    if reentries[i] == 0
        reentries[i] = 1e-5
        years_out_before_reentry[i] = 1
    end

    # Frac Immediate Re-Entry
    exits_all_noEndYears[i] = sum(exits_all[2:11, :, :])
    immediate_reentries[i] = sum(Int.((exits_all[2:11, :, :] .== 1) .& (entries_all[3:12, :, :] .== 1)))

    # If there are no exits, change it such that the data will be very far off the moment
    if exits_all_noEndYears[i] == 0
        exits_all_noEndYears[i] = 1e-5
        immediate_reentries[i] = 1
    end
end
avg_years_out_delta = years_out_before_reentry./reentries
plot(avg_years_out_delta)
frac_immediate_reentry_delta = immediate_reentries./exits_all_noEndYears
plot(frac_immediate_reentry_delta)

row = 1
for k = 1:prim.n_sims
    for j = 1:prim.n_firms
        for i = 1:12
            # Panel will be: Firm Year Export Capital Sales Export Sales
            panel[row,1] = parse(Float32, string(k,".",j,1))
            panel[row,2] = i
            panel[row,3] = firm_export_choices[prim.n_periods_experiment-12+i,j,k]
            panel[row,4] = firm_capital_choices[prim.n_periods_experiment-12+i,j,k]
            panel[row,5] = firm_sales_all[prim.n_periods_experiment-12+i,j,k]
            panel[row,6] = firms_export_sales[prim.n_periods_experiment-12+i,j,k]
            panel[row,7] = productivities[prim.n_periods_experiment-12+i,j,k]
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

###
# Create %Δ Plots for Canonical
###
lchange_sales_Base = zeros(prim.n_periods_experiment-1, prim.n_firms, prim.n_sims)
lchange_sales_Base = (firm_sales_all[2:prim.n_periods_experiment,:,:].-firm_sales_all[1:prim.n_periods_experiment-1,:,:])./firm_sales_all[1:prim.n_periods_experiment-1,:,:]
avg_lchange_sales_Base = mean(lchange_sales_Base, dims = [2,3])
plot(avg_lchange_sales_Base[2:prim.n_periods_experiment-1])

total_export_sales_Base = sum(firms_export_sales, dims = 2)
plot(mean(total_export_sales_Base,dims=3)[:,1,1])

total_sales_Base = sum(firm_sales_all, dims = 2)
plot(mean(total_sales_Base,dims=3)[:,1,1])

lchange_export_sales_Base = zeros(prim.n_periods_experiment-1, prim.n_firms, prim.n_sims)
lchange_export_sales_Base = (firms_export_sales[2:prim.n_periods_experiment,:,:].-firms_export_sales[1:prim.n_periods_experiment-1,:,:])./firms_export_sales[1:prim.n_periods_experiment-1,:,:]
avg_lchange_export_sales_Base = mean(lchange_export_sales_Base, dims = [2,3])
plot(avg_lchange_export_sales_Base[2:prim.n_periods_experiment-1])

firms_export_decisions_no_zero = hcat(firm_export_choices, ones(prim.n_periods_experiment, 1, prim.n_sims), zeros(prim.n_periods_experiment, 1, prim.n_sims)) 
firms_export_change = zeros(prim.n_periods_experiment, prim.n_firms, prim.n_sims)
for i = 2:prim.n_periods_experiment
    for j = 1:prim.n_firms
        firms_export_change[i,j,:] = firm_export_choices[i,j,:] .- firm_export_choices[i-1,j,:]
    end
end
firms_export_change_no_zero = hcat(firms_export_change, ones(prim.n_periods_experiment, 1, prim.n_sims), -1 .* ones(prim.n_periods_experiment, 1, prim.n_sims)) 
starter_rate = zeros(prim.n_periods_experiment-1,prim.n_sims)
stopper_rate = zeros(prim.n_periods_experiment-1,prim.n_sims)

exporter_i = countmap(firms_export_change_no_zero[1,:,:])
exporter_i_minus1 = countmap(firms_export_change_no_zero[1,:,:])
for i = 2:prim.n_periods_experiment-1
    for k = 1:prim.n_sims
        exporter_i_minus1 = countmap(firms_export_decisions_no_zero[i-1,:,k])
        exporter_i = countmap(firms_export_change_no_zero[i,:,k])
        
        starter_rate[i-1,k] = exporter_i[1]/exporter_i_minus1[0]
        stopper_rate[i-1,k] = exporter_i[-1]/exporter_i_minus1[1]
    end
end
avg_starter_rate_Base = mean(starter_rate, dims=2)
plot(avg_starter_rate_Base)
avg_stopper_rate_Base = mean(stopper_rate,dims=2)
plot(avg_stopper_rate_Base)

## Coef of Variation Over Time
coef_variation = zeros(prim.n_periods_experiment-20, prim.n_sims)
for j = 1:prim.n_sims
    for i = 1:prim.n_periods_experiment-20
        coef_variation[i,j] = std(log.(firm_sales_domestic[i:i+20,:,j]))/mean(log.(firm_sales_domestic[i:i+20,:,j]))
    end
end
coef_variation_Base = mean(coef_variation,dims=2)
plot(coef_variation_Base)

## β-Coef Over Time
β_sales_Base = zeros(prim.n_periods_experiment-20, 1)
for i = 2:prim.n_periods_experiment-20
    for s = 1:1
        sales_today = firm_sales_domestic[i+1,:,s]
        sales_yesterday = firm_sales_domestic[i,:,s]
        year = []
        firm_num = []
        for i = i+2:i+19
            sales_today = [sales_today; firm_sales_domestic[i,:,s]]
            sales_yesterday = [sales_yesterday; firm_sales_domestic[i-1,:,s]]
        end

        for i = 2:20
            for j = 1:prim.n_firms
                year = [year; string("Y",i)]
                firm_num = [firm_num; string("F",j)]
            end
        end

        sales_today_vec = copyto!(Vector{Float64}(undef,length(sales_today)),sales_today)
        sales_yesterday_vec = copyto!(Vector{Float64}(undef,length(sales_yesterday)),sales_yesterday)

        reg_DF = DataFrame(s_today = log.(sales_today_vec), s_yesterday = log.(sales_yesterday_vec), y = year, f = firm_num)
        reg_results = reg(reg_DF, @formula(s_today ~ s_yesterday + fe(y) + fe(f)))

        β_moment = coef(reg_results)
        β_sales_Base[i,s] = β_moment[1]
    end
end
plot(β_sales_Base)

## Export-Sales Ratio Over Time
export_sales_ratio_Base = zeros(prim.n_periods_experiment-20)
for i = 1:prim.n_periods_experiment-20
    export_sales_ratio_Base[i] = sum(firms_export_sales[i:i+19,:,:]./firm_sales_all[i:i+19,:,:])/sum(firm_export_choices[i:i+19,:,:])
end
plot(export_sales_ratio_Base)

## Average Years Out and Frac Immedate Re-Entry Over Time
avg_years_out_Base = zeros(prim.n_periods_experiment-12)
frac_immediate_reentry_Base = zeros(prim.n_periods_experiment-12)
reentries = zeros(prim.n_periods_experiment-12)
years_out_before_reentry = zeros(prim.n_periods_experiment-12)
exits_all_noEndYears = zeros(prim.n_periods_experiment-12)
immediate_reentries = zeros(prim.n_periods_experiment-12)

for i=1:prim.n_periods_experiment-12
    
    exits_all = Int.(firms_export_change[i:i+11,:,:] .== -1)
    exits_cumulative = Int.(firms_export_change[i:i+11,:,:] .== -1)
    final_exit_no_reentry = ones(Int, 12, prim.n_firms, prim.n_sims)
    entries_all = Int.(firms_export_change[i:i+11,:,:] .== 1)

    for j in 2:12
        exits_cumulative[j,:,:] .+= exits_cumulative[j-1,:,:]
    end

    # Set entries where the firm is exporting equal to zero, keeping strings of numbers with length equal to years out of each exit
    exits_cumulative = exits_cumulative .* (ones(Int, 12, prim.n_firms, prim.n_sims) .- firm_export_choices[i:i+11,:,:])
    # Test if that exit is the final exit and has no re-entry tied to it
    for j in 1:11
        final_exit_no_reentry[j,:,:] = Int.((exits_cumulative[j,:,:] .== exits_cumulative[12,:,:]))
    end
    # Only keep exits that result in re-entries
    exits_cumulative = exits_cumulative .* (1 .- final_exit_no_reentry)
    # Loop through the number of exits and sum up the years out that result in re-entry
    max_exits = maximum(exits_cumulative)
    exits_count_map = countmap(exits_cumulative)
    years_out_before_reentry[i] = 0
    for j in 1:Int(max_exits)
        # Test if that number of exits actually occurs ever
        try
            exits_count_map[j]
        catch
            # if not set it to 0
            exits_count_map[j] = 0
        end
        years_out_before_reentry[i] += exits_count_map[j]
    end
    # Get maximum exits for each firm-simulation, which is number of re-entries
    reentries[i] = sum(maximum(exits_cumulative, dims = 1))
    # If there are no re-entries make them very small so it will be very far off the moment
    if reentries[i] == 0
        reentries[i] = 1e-5
        years_out_before_reentry[i] = 1
    end

    # Frac Immediate Re-Entry
    exits_all_noEndYears[i] = sum(exits_all[2:11, :, :])
    immediate_reentries[i] = sum(Int.((exits_all[2:11, :, :] .== 1) .& (entries_all[3:12, :, :] .== 1)))

    # If there are no exits, change it such that the data will be very far off the moment
    if exits_all_noEndYears[i] == 0
        exits_all_noEndYears[i] = 1e-5
        immediate_reentries[i] = 1
    end
end
avg_years_out_Base = years_out_before_reentry./reentries
plot(avg_years_out_Base)
frac_immediate_reentry_Base = immediate_reentries./exits_all_noEndYears
plot(frac_immediate_reentry_Base)

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