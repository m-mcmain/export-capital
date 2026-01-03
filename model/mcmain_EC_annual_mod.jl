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
# x0_first3 = [1.7353879062864415, 0.47, 0.14254513790307305, 0.71, 0.16898585564776006]
# # [1.3056182993255299, 0.41040338483719635, 0.14315240437674726, 0.16457325094720668, 0.7390948399070301] first best for w = 0.25
# #[2.031499422574203, 0.4691582010613471, 0.195338159196638] #, 0.825, 0.07047818746078156] # error of 0.002
# opt_res_first3 = optimize(MSM_func_first3, x0_first3)
# minimizers_first3 = Optim.minimizer(opt_res_first3)
# minimizers_first3

# just for Optim, comment out below:


##############################################################
#####                     Optim Delta                     ####
##############################################################
prim, res = Initialize(3)
# [0.04672832247947507, 0.00016618510618323404, 5.208544661781648, 0.4768979634848714],,,, 0.63 & 0.3
x0_delta_first3 = [0.02659622172414043, 5.020533056952626, 0.4562355246397134, 0.14315240437674726, 0.7099550149528832, 0.17358427873116328] 
# [0.0005199261953161846, 0.0009185731464601314, 0.02659622172414043, 5.020533056952626, 0.4562355246397134, 0.14315240437674726, 0.7099550149528832, 0.17358427873116328] final one i saw, moving to just \delta
# [0.006126027377287263, 0.00879566685351088, 0.013400275982553203, 5.969505448035115, 0.4400303853182756] try 3, similar? Better than 1 by a bit
# [0.0008894943383066979, 0.0202801249481327, 0.00947766961379277, 6.037713174594099, 0.40360335262628455] try 2, similar results?
# [0.0005199261953161846, 0.0009185731464601314, 0.02659622172414043, 5.020533056952626, 0.4562355246397134] try 1
# [0.024369523724071876, 0.047216604867249844, 5.742159276459619, 0.45247753941501945]
#[0.03481335347534838, 5.748033462346467, 0.440448204725421] # [0.04, 3.249101472920794, 0.5]#[0.1, 1.5, 0.2]
# [0.05583548816343635, 4.382436294248027, 0.45159997957845405] different best, 0.69 & 0.275
#[0.07687879267482975, 3.8541397148948935, 0.45855335601134445] different best, fixed re-entry %s
#[0.007744847458919135, 9.874385277937614, 0.43464896718069285] first best?
#[0.018464576233428683, 5.108017946656351, 0.40360485919683664] second best?
#[0.030260316445346234, 0.6156891515679773, 8.28542705594612] # Fixed Cost based on Last Export

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