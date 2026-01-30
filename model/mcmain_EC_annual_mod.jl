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
end

using Parameters, Plots, Optim, Distributions, SharedArrays, Distributed, Random, JLD2, Statistics, StatsBase, GLM, DataFrames, OrderedCollections, LinearAlgebra, FixedEffectModels, QuantEcon

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
prim, res = Initialize(3)
# [0.04672832247947507, 0.00016618510618323404, 5.208544661781648, 0.4768979634848714],,,, 0.63 & 0.3
x0_delta_first3 = [0.008489760052962878, 0.02553618111567095, 8.273244642640746, 0.6813261628529305, 0.23764509440810713]#, 0.7099550149528832, 0.172] 
# [0.009967668768297765, 0.06310863813742076, 6.693685173035949, 0.6983557136908796, 0.23764510559905438]
# Good on all but E/S Ratio 0.01 higher and Years Out 1.7 vs. 1.95, unmatched is meh
# [0.008703447586350862, 0.10344459793805771, 6.813178391349231, 0.735309605793568, 0.25048947468382765, 0.7114245033293098, 0.17957495593216935] # gone 5 years out of market
# [0.010067553953214091, 0.2760026539936621, 5.383209003828556, 0.70598793627047, 0.23764510550576487] gone 3 years out of market

# [0.005115859873079168, 23.775353858407993, 0.49993628089731706, 0.23764510272853265] Error: 0.08797
# [0.009950719651712845, 0.050624862896243086, 6.850584868430243, 0.6124504044678348, 0.23764510620543247] 
opt_res_delta_first3 = optimize(MSM_delta_func_first3, x0_delta_first3)
minimizers_delta_first3 = Optim.minimizer(opt_res_delta_first3)
minimizers_delta_first3



#=
#x0_last2 = [0.8027153603452549, 0.11175754485806927]
#lower = [0.01, 0.01]
#upper = [0.99, Inf]
#inner_optimizer = GradientDescent()
#opt_res_last2 = optimize(MSM_func_last2, lower, upper, x0_last2, Fminbox(inner_optimizer))
#minimizers_last2 = Optim.minimizer(opt_res_last2)
#print(minimizers_last2)

##############################################################
#####              Test with close Params                 ####
##############################################################
# Make sure to put in correct values of \delta, FC_0, and FC_1
# [0.1412008592544308, 0.4595257572507812, 6.18159395178584]
# other close is [0.085, 0.535, 6.55]
prim, res = Initialize(3)
Solve_model(prim, res)

@unpack val_func, ex_func = res
@unpack ϵ_grid = prim
Plots.plot(ϵ_grid, ex_func[3, :, :], title="Export Function, Q=1",
label=["0" "d^11" "d^10" "d^9" "d^8" "d^7" "d^6" "d^5" "d^4" "d^3" "d^2" "d" "1"], linewidth=3, xlabel="ϵ")
savefig("Epsilon_Cutoffs.png")

firms_export_decisions, firms_labor_decisions, firms_capital_decisions, firms_sales_domestic, firms_sales = data_sim_delta_nsims(prim, res)
prop_one_entry_mean, prop_two_entry_mean, prop_three_entry_mean, prop_four_entry_mean, prop_five_entry_mean = reentry_calcs(firms_export_decisions, firms_sales_domestic, firms_sales)
bar([1,2,3,4,5],[prop_one_entry_mean, prop_two_entry_mean, prop_three_entry_mean, prop_four_entry_mean, prop_five_entry_mean], legend = false, dpi=300)
xlabel!("Number of Entries")
ylabel!("Percent")
title!(" ")
annotate!((1.5,prop_one_entry_mean,(prop_one_entry_mean,:bottom,11)))
annotate!((2.5,prop_two_entry_mean,(prop_two_entry_mean,:bottom,11)))
annotate!((3.5,prop_three_entry_mean,(prop_three_entry_mean,:bottom,11)))
annotate!((4.5,prop_four_entry_mean,(prop_four_entry_mean,:bottom,11)))
annotate!((5.5,prop_four_entry_mean,(prop_five_entry_mean,:bottom,11)))
savefig("Figure1_Calibrated_Model_trim.png")
=#