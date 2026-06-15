@with_kw struct Primitives
    σ_q::Float64 = 0.07 # s.e. for RER shock
    ρ_q::Float64 = 0.835 # AR coef for RER shock
    R::Float64 = 1/(1+0.124) # Discount Rate (annualized is 12.4%)
    r::Float64 = 0.124 # Interest Rate
    w::Float64 = 0.225 # Wage guess?
    α_n::Float64 = 0.505 # Labor share of income
    α_k::Float64 = 0.495 # Plant level returns to scale
    θ::Float64 = 5.0 # Elasticity of Substitution
    C::Float64 = 1.0 # Domestic Aggregate Consumption
    γ_0::Float64 = 0.258 # entry level of demand
    γ_1::Float64 = 0.024 # demand gain from additional year in market
        
    Q_min::Float64 = 0.8 # RER lower bound
    Q_max::Float64 = 1.2 # RER upper bound
    nQ::Int64 = 5 #number of capital grid points
    Q_grid::Array{Float64,1} = collect(range(Q_min, length = nQ, stop = Q_max))  # RER grid
    
    ϵ_min::Float64 = 0.7 # Productivity Lower Bound
    ϵ_max::Float64 = 1.4 # Productivity Upper Bound
    nϵ::Int64 = 101 #Number of Productivity states 
    ϵ_grid::Array{Float64,1} = collect(range(ϵ_min, length = nϵ, stop = ϵ_max)) # Productivity Grid

    n_periods::Int64 = 112 # Number of periods
    n_periods_experiment::Int64 = 112 # Number of periods
    n_firms::Int64 = 2826 # Number of Firms
    n_sims::Int64 = 1000 # Number of simulations

    true_starter::Float64 = 0.10091176 # True starter rate
    true_stopper::Float64 = 0.11991836 # True stopper rate
    true_ave_es_ratio::Float64 = 0.18197473 # True average exports to sales ratio
    true_enter_twice_perc::Float64 = 0.2548 # True firms who enter twice
    true_enter_once_perc::Float64 = 0.6042 # True firms who enter once
    true_enter_thrice_perc::Float64 = 0.1171 # True firms who enter thrice
    true_enter_four_perc::Float64 = 0.022 # True firms who enter four times
    true_enter_five_perc::Float64 = 0.0018 # True firms who enter five times
    true_coef_var::Float64 = 0.09592415 # True coefficient of variation
    true_a_exp_growth::Float64 = 0.40958202 # True regression coefficient on log sales growth
    true_avg_time_out_reentrant::Float64 = 1.9461679 # True average time out of the export market for re-entrants
    true_perc_reenter_immediate::Float64 = 0.4069149 # True percent of firms that reenter after 1 year out

end

mutable struct Results
    n_func::Array{Float64, 3} # Optimal Workers
    k_func::Array{Float64, 3} # Optimal Capital
    ex_func::Array{Float64, 3} # Optimal Export Decision
    ex_cap::Array{Float64, 2} # Export Capital
    ϵ::Array{Float64, 3} # Productivity for the firm
    Q::Array{Float64, 2} # Real Exchange Rate
    ϵ_experiment::Array{Float64, 2} # Productivity for the firm
    Q_experiment::Array{Float64, 1} # Real Exchange Rate
    val_func::Array{Float64, 3} # Value function
    τ::Float64
    # Estimated Parameters:
    C_star::Float64 # Foreign demand scale guess
    σ_e::Float64 # s.e. for productivity shock guess
    ρ_e::Float64 # AR coef for productivity shock guess
    FC_0::Float64 # Sunk Cost based on Last Export
    FC_1::Float64 # Fixed Cost reentry based on Last Export
    # Multiplied by median sales of non_exporters w/ these parameters (18)
    δ::Float64 # Export Knowledge Depreciation
    α_d::Float64 # Decay function intercept
    β_d::Float64 # Decay function coefficient
    β_sq_d::Float64 # Decay function squared coefficient
    n_prev_ex::Int64 # Previous export states
    prev_ex_grid::Array{Float64,1} # Previous export grid
    tauchen_trans_Q::Array{Float64,2} # Tauchen's Method Transition Probs for Q
    tauchen_trans_e::Array{Float64,2} # Tauchen's Method Transition Probs for ϵ
    # Res for Subsidy Experiment
    val_func_subsidy::Array{Float64, 3} # Value function with one time subsidy
    n_func_subsidy::Array{Float64, 3} # Optimal Workers
    k_func_subsidy::Array{Float64, 3} # Optimal Capital
    ex_func_subsidy::Array{Float64, 3} # Optimal Export Decision
end

function decay_grid(intercept::Float64, coefficient::Float64, coef_sq::Float64, n_grid::Int64, type::Int64)
    # Type = 1 means Quadratic. Exponential, otherwise. 

    if type == 1
        grid_d = [1.0]
        for i =0:n_grid-1
            d_left = min(max(1-(intercept + coefficient*i+coef_sq*i^2),0),1)
            append!(grid_d, d_left)
        end
    else
        grid_d = [1.0]
        for i = 0:n_grid-1
            d_left = min(max(1-(intercept+coefficient*i^2),0),1)
            append!(grid_d, d_left)
        end
    end

    append!(grid_d, 0)
    return reverse(grid_d)
end

function Initialize(base::Int64)
    prim = Primitives() #initialize primtiives
    ex_cap = Array{Float64}(zeros(prim.n_periods, prim.n_firms)) # initial export capital
    ϵ = Array{Float64}(ones(prim.n_periods, prim.n_firms, prim.n_sims)) # initial ϵ
    Q = Array{Float64}(ones(prim.n_periods, prim.n_sims)) # initial Q
    ϵ_experiment = Array{Float64}(ones(prim.n_periods_experiment, prim.n_firms)) # initial ϵ
    Q_experiment = Array{Float64}(ones(prim.n_periods_experiment)) # initial Q
    τ = 0.0 # initial τ (no tariffs)
    if base == 1
        C_star = 0.14822662063483324
        ρ_e =  0.5704358764309491
        σ_e = 0.05389818114125237
        FC_0 = 1.6741901913074753
        FC_1 = 0.4352053721156689
        δ = 1
        α_d = 1.0 # Export Capital Decay Intercept Guess 
        β_d = 0.0 # Export Capital Decay coefficient Guess       
        β_sq_d = 0.0 # Export Capital Decay squared coefficient Guess  
    elseif base == 2 # Perfectly Identified Sunk Cost Model 
        C_star = 0.13017264011231047
        ρ_e =  0.9999843253276733
        σ_e = 0.008942722798834544
        FC_0 = 0.08176178216576599
        FC_1 = 0.8307762231237128
        δ = 1
        α_d = 1.0 # Export Capital Decay Intercept Guess 
        β_d = 0.0 # Export Capital Decay coefficient Guess       
        β_sq_d = 0.0 # Export Capital Decay squared coefficient Guess  
    elseif base == 3
        C_star = 0.16333309353877032
        ρ_e =  0.5717084452766545
        σ_e = 0.05447497950369329
        FC_0 = 3.0133789036830056 # Fixed Cost based on Last Export
        FC_1 = 0.5286860069055453 # Fixed Cost of reentry
        δ = 0.03379118541575009 # Export Knowledge Deprication guess 
        α_d = 0.03379118541575009 # Export Capital Decay Intercept Guess 
        β_d = 0.003269279322491091 # Export Capital Decay coefficient Guess       
        β_sq_d = 0.0 # Export Capital Decay squared coefficient Guess    
    end

    if δ == 1
        n_prev_ex = 2 # Number of previous export states
        prev_ex_grid = [0, 1]
    else
        prev_ex_grid = unique(decay_grid(α_d, β_d, β_sq_d, 22-2, 2))
        n_prev_ex = length(prev_ex_grid)
    end
    n_func = Array{Float64}(zeros(prim.nQ, prim.nϵ, n_prev_ex)) # initial worker function guess
    k_func = Array{Float64}(zeros(prim.nQ, prim.nϵ, n_prev_ex)) # initial capital function guess
    ex_func = Array{Float64}(zeros(prim.nQ, prim.nϵ, n_prev_ex)) # initial export function guess
    val_func = Array{Float64}(zeros(prim.nQ, prim.nϵ, n_prev_ex)) # initial value function guess

    # Res for Subsidy Experiment
    n_func_subsidy = Array{Float64}(zeros(prim.nQ, prim.nϵ, n_prev_ex)) # initial worker function guess
    k_func_subsidy = Array{Float64}(zeros(prim.nQ, prim.nϵ, n_prev_ex)) # initial capital function guess
    ex_func_subsidy = Array{Float64}(zeros(prim.nQ, prim.nϵ, n_prev_ex)) # initial export function guess
    val_func_subsidy = Array{Float64}(zeros(prim.nQ, prim.nϵ, n_prev_ex)) # initial value function guess

    tauchen_res_Q = tauchen(prim.nQ, prim.ρ_q, prim.σ_q)
    tauchen_trans_Q = tauchen_res_Q.p

    tauchen_res_e = tauchen(prim.nϵ, ρ_e, σ_e)
    tauchen_trans_e = tauchen_res_e.p

    res = Results(n_func, k_func, ex_func, ex_cap, ϵ, Q, ϵ_experiment, Q_experiment, val_func, τ, C_star, σ_e, ρ_e, FC_0, FC_1, δ, α_d, β_d, β_sq_d, n_prev_ex, prev_ex_grid, tauchen_trans_Q, tauchen_trans_e, val_func_subsidy, ex_func_subsidy, n_func_subsidy, k_func_subsidy) #initialize results struct
    prim, res #return deliverables
end

#Bellman Operator
function Bellman(prim::Primitives,res::Results)
    @unpack val_func, ex_func, n_func, k_func, τ, C_star, σ_e, ρ_e, FC_0, FC_1, prev_ex_grid, n_prev_ex, tauchen_trans_Q, tauchen_trans_e = res #unpack value function
    @unpack Q_grid, Q_max, ϵ_grid, ϵ_max, R, r, w, σ_q, nQ, Q_grid, nϵ, ϵ_grid, ρ_q, α_n, θ, C, γ_0, γ_1 = prim #unpack model primitives
    v_next = zeros(nQ,nϵ,n_prev_ex) #next guess of value function to fill

    for Q_index = 1:nQ   
        Q = Q_grid[Q_index] #value of Q
        for ϵ_index = 1:nϵ
            ϵ = ϵ_grid[ϵ_index] #value of ϵ
            for prev_ex_index in 1:n_prev_ex
                prev_ex = prev_ex_grid[prev_ex_index]
                candidate_max = -1e10 #bad candidate max
                
                n_prev = 0
                k_prev = 0
                Π_prev = 0

                for ex_decision in 0:1 #loop over possible export decision
                    Π(x) = -1*((1+ex_decision*(Q*(1-res.τ))^θ*C_star/C)^(1/θ)*C^(1/θ)*ϵ*x[1]^(α_n*(θ-1)/θ)*x[2]^((1-α_n)*(θ-1)/θ)-w*x[1]-r*x[2])
                    x0 = [0.1, 0.1]
                    lower = [0.01, 0.01]
                    upper = [Inf, Inf]
                    inner_optimizer = GradientDescent()
                    optim_Π = optimize(Π, lower, upper, x0, Fminbox(inner_optimizer))

                    n_opt = Optim.minimizer(optim_Π)[1]
                    k_opt = Optim.minimizer(optim_Π)[2]
                    Π = -1*Optim.minimum(optim_Π)

                    # Should this be different ? Adjust based on export choice?
                    continuation_val = 0
                    for ϵp_index in 1:nϵ
                        for Qp_index in 1:nQ
                            continuation_val += R*tauchen_trans_e[ϵ_index,ϵp_index]*tauchen_trans_Q[Q_index,Qp_index]*val_func[Qp_index,ϵp_index,prev_ex_index]
                        end
                    end
                    total_Π = Π + continuation_val
                    if ex_decision == 0
                        Π_prev = total_Π
                        n_prev = n_opt
                        k_prev = k_opt
                    else
                        if prev_ex_index == n_prev_ex
                            total_Π += -1*FC_1
                        else
                            total_Π += -1*((1-prev_ex)*FC_0+FC_1)
                        end

                        if total_Π > Π_prev
                            ex_func[Q_index, ϵ_index, prev_ex_index] = 1
                            n_func[Q_index, ϵ_index, prev_ex_index] = n_opt
                            k_func[Q_index, ϵ_index, prev_ex_index] = k_opt
                            if total_Π > candidate_max
                                candidate_max = total_Π
                            end
                        else
                            ex_func[Q_index, ϵ_index, prev_ex_index] = 0
                            n_func[Q_index, ϵ_index, prev_ex_index] = n_prev
                            k_func[Q_index, ϵ_index, prev_ex_index] = k_prev
                            if Π_prev > candidate_max
                                candidate_max = Π_prev
                            end
                        end
                    end

                end
                v_next[Q_index, ϵ_index, prev_ex_index] = candidate_max #update value function
            end
        end
    end
    v_next #return next guess of value function
end

#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-3, err::Float64 = 1000.0)
    n = 0 #counter

    while err>tol #begin iteration
        v_next = Bellman(prim, res) #spit out new vectors
        err = abs.(maximum(v_next.-res.val_func))/abs(v_next[prim.nQ, prim.nϵ, 1]) #reset error level
        res.val_func = v_next #update value function
        #println("Error:", err)
        n+=1
    end
    println("Value function converged in ", n, " iterations.")
end

#solve the model
function Solve_model(prim::Primitives, res::Results)
    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
end

function MSM_func_first3(x)
    prim, res = Initialize(1) #initialize primitive and results structs

    res.FC_0 = x[1]
    res.FC_1 = x[2]
    res.C_star = x[3]
    res.ρ_e = x[4]
    res.σ_e = x[5]
    

    @time V_iterate(prim, res)

    firms_export_decisions, firms_labor_decisions, firms_capital_decisions, firms_sales_domestic, firms_sales = data_sim_delta_nsims(prim, res)

    val_annual_firms_export_decisions = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_export_decisions = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_export_change = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales_domestic = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales_foreign = zeros(12, prim.n_firms, prim.n_sims)
    annual_total_sales_foreign = zeros(12, prim.n_sims)
    annual_total_sales = zeros(12)

    for j = 1:12
        val_annual_firms_export_decisions[j,:,:] = val_annual_firms_export_decisions[j,:,:] .+ firms_export_decisions[j+100,:,:]
        annual_firms_sales_domestic[j,:,:] = annual_firms_sales_domestic[j,:,:] .+ firms_sales_domestic[j+100,:,:]
        annual_firms_sales[j,:,:] = annual_firms_sales[j,:,:] .+ firms_sales[j+100,:,:]
        annual_firms_sales_foreign[j,:,:] = annual_firms_sales_foreign[j,:,:] .+ firms_sales[j+100,:,:] .- firms_sales_domestic[j+100,:,:]
    end
       
    annual_firms_export_decisions = val_annual_firms_export_decisions .> 0
    mean_annual_firms_sales_foreign = mean(annual_firms_sales_foreign, dims=3)
    mean_annual_firms_sales = mean(annual_firms_sales, dims=3)

    for i = 1:12
        for j = 1:prim.n_firms
            if i == 1
                annual_total_sales_foreign[i] = annual_total_sales_foreign[i] + mean_annual_firms_sales_foreign[i,j]
                annual_total_sales[i] = annual_total_sales[i] + mean_annual_firms_sales[i,j]
            else
                annual_total_sales_foreign[i] = annual_total_sales_foreign[i] + mean_annual_firms_sales_foreign[i,j] 
                annual_total_sales[i] = annual_total_sales[i] + mean_annual_firms_sales[i,j]
                annual_firms_export_change[i,j,:] = annual_firms_export_decisions[i,j,:] .- annual_firms_export_decisions[i-1,j,:]
            end
        end
    end
    annual_starter_rate = zeros(11, prim.n_sims)
    annual_stopper_rate = zeros(11, prim.n_sims)
    # Add 1 and 0 to export decisions once to make sure neither are ever 0
    annual_firms_export_decisions = hcat(annual_firms_export_decisions, ones(12, 1, prim.n_sims), zeros(12, 1, prim.n_sims))
    # Add 1 and -1 to export changes once to make sure neither are ever 0
    annual_firms_export_change = hcat(annual_firms_export_change, ones(12, 1, prim.n_sims), -1 .* ones(12, 1, prim.n_sims)) 
    exporter_i = countmap(annual_firms_export_change[1,:,:])
    exporter_i_minus1 = countmap(annual_firms_export_decisions[1,:,:])

    for i = 2:12
        for k = 1:prim.n_sims
            exporter_i_minus1 = countmap(annual_firms_export_decisions[i-1,:,k])
            exporter_i = countmap(annual_firms_export_change[i,:,k])
            
            annual_starter_rate[i-1,k] = exporter_i[1]/exporter_i_minus1[0]
            annual_stopper_rate[i-1,k] = exporter_i[-1]/exporter_i_minus1[1]
        end
    end

    sales_today = mean(annual_firms_sales[2,:,:], dims=2)
    sales_yesterday = mean(annual_firms_sales[1,:,:], dims=2)
    year = []
    firm_num = []
    for i = 3:12
        sales_today = [sales_today; mean(annual_firms_sales[i,:,:],dims=2)]
        sales_yesterday = [sales_yesterday; mean(annual_firms_sales[i-1,:,:],dims=2)]
    end

    for i = 2:12
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

    output = zeros(5)

    output[1] = mean(annual_starter_rate)
    output[2] = mean(annual_stopper_rate)
    output[3] = mean(res.C_star ./ (res.Q[101:112,17].^(-1*prim.θ)))
    output[4] = std(log.(annual_firms_sales_domestic[:,:,20]))/mean(log.(annual_firms_sales_domestic[:,:,20]))
    output[5] = β_moment[1]

    error = abs.(output[1]-prim.true_starter)+abs.(output[2]-prim.true_stopper)+abs.(output[4]-prim.true_coef_var)+abs.(output[5]-prim.true_a_exp_growth)+abs.(output[3]-prim.true_ave_es_ratio)
    #
    print(error)
    print("\n")
    print(x)
    print("\n")
    print(output)
    print("\n")
    print([prim.true_starter, prim.true_stopper, prim.true_ave_es_ratio, prim.true_coef_var, prim.true_a_exp_growth]) 
    print("\n")
    print(median(annual_firms_sales))
    print("\n")
    print(median(annual_firms_sales_foreign))
    print("\n")
    print(median(firms_labor_decisions[101:112,:,:]))
    print("\n")
    print(median(firms_capital_decisions[101:112,:,:]))
    print("\n")

    return error 

end

function MSM_delta_func_first3(x)
    print(x)
    print("\n")
    model = 3
    prim, res = Initialize(model) #initialize primitive and results structs

    if model == 3
        res.α_d = x[1]
        res.β_d = x[2]

        # if x[6] >= 1
        #     x[6] = 0.99
        # end

        # if x[7] < 0
        #     x[7] = 0.01
        # end

        # res.C_star = x[5]
        # res.ρ_e = x[6]
        # res.σ_e = x[7]

    else
        if x[1] >= 1
            x[1] = 0.99
        end

        if x[2] < 0
            x[2] = 0.01
        end

        res.C_star = x[5]
        res.ρ_e = x[1]
        res.σ_e = x[2]
        print([res.ρ_e res.σ_e])
    end

    # FC_0 is x[3] no matter what
    if x[3] < 0
        x[3] = 0
    end
    res.FC_0 = x[3]

    # # FC_1 is x[4] no matter what
    if x[4] < 0
        x[4] = 0
    end
    res.FC_1 = x[4]

    # Update the ϵ process
    # tauchen_res_e = tauchen(prim.nϵ, res.ρ_e, res.σ_e)
    # res.tauchen_trans_e = tauchen_res_e.p

    # Update the export decay grid
    if res.δ < 1
        res.prev_ex_grid = unique(decay_grid(x[1], x[2], x[3], 24-2, 2))
        res.n_prev_ex = length(res.prev_ex_grid)
    else
        res.prev_ex_grid = [0, 1]
        res.n_prev_ex = length(res.prev_ex_grid)
    end

    print(res.prev_ex_grid)

    # Initialize functions
    res.n_func = Array{Float64}(zeros(prim.nQ, prim.nϵ, res.n_prev_ex)) # initial worker function guess
    res.k_func = Array{Float64}(zeros(prim.nQ, prim.nϵ, res.n_prev_ex)) # initial capital function guess
    res.ex_func = Array{Float64}(zeros(prim.nQ, prim.nϵ, res.n_prev_ex)) # initial export function guess
    res.val_func = Array{Float64}(zeros(prim.nQ, prim.nϵ, res.n_prev_ex)) # initial value function guess
    @time V_iterate(prim, res)

    firms_export_decisions, firms_labor_decisions, firms_capital_decisions, firms_sales_domestic, firms_sales = data_sim_delta_nsims(prim, res)
    median(firms_labor_decisions[101:112,:,:])
    
    val_annual_firms_export_decisions = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_export_decisions = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_export_change = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales_domestic = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales_foreign = zeros(12, prim.n_firms, prim.n_sims)
    annual_total_sales_foreign = zeros(12, prim.n_sims)
    annual_total_sales = zeros(12)
    prop_one_entry = zeros(prim.n_sims)
    prop_two_entry = zeros(prim.n_sims)
    prop_three_entry = zeros(prim.n_sims)
    prop_four_entry = zeros(prim.n_sims)
    prop_five_entry = zeros(prim.n_sims)
    prop_six_entry = zeros(prim.n_sims)
    prop_seven_entry = zeros(prim.n_sims)
    prop_no_entry_all = zeros(prim.n_sims)
    prop_one_entry_all = zeros(prim.n_sims)
    prop_two_entry_all = zeros(prim.n_sims)
    prop_three_entry_all = zeros(prim.n_sims)
    prop_four_entry_all = zeros(prim.n_sims)
    avg_time_out_reentrant = zeros(prim.n_sims)

    for j = 1:12
        val_annual_firms_export_decisions[j,:,:] = val_annual_firms_export_decisions[j,:,:] .+ firms_export_decisions[j+100,:,:]
        annual_firms_sales_domestic[j,:,:] = annual_firms_sales_domestic[j,:,:] .+ firms_sales_domestic[j+100,:,:]
        annual_firms_sales[j,:,:] = annual_firms_sales[j,:,:] .+ firms_sales[j+100,:,:]
        annual_firms_sales_foreign[j,:,:] = annual_firms_sales_foreign[j,:,:] .+ firms_sales[j+100,:,:] .- firms_sales_domestic[j+100,:,:]
    end
       
    annual_firms_export_decisions = val_annual_firms_export_decisions .> 0
    annual_firms_export_change[1,:,:] = copy(annual_firms_export_decisions[1,:,:])
    mean_annual_firms_sales_foreign = mean(annual_firms_sales_foreign, dims=3)
    mean_annual_firms_sales = mean(annual_firms_sales, dims=3)

    for i = 1:12
        for j = 1:prim.n_firms
            if i == 1
                annual_total_sales_foreign[i] = annual_total_sales_foreign[i] + mean_annual_firms_sales_foreign[i,j]
                annual_total_sales[i] = annual_total_sales[i] + mean_annual_firms_sales[i,j]
            else
                annual_total_sales_foreign[i] = annual_total_sales_foreign[i] + mean_annual_firms_sales_foreign[i,j] 
                annual_total_sales[i] = annual_total_sales[i] + mean_annual_firms_sales[i,j]
                annual_firms_export_change[i,j,:] = annual_firms_export_decisions[i,j,:] .- annual_firms_export_decisions[i-1,j,:]
            end
        end
    end

    ########################################
    # Calculate starter and stopper rates  #
    ########################################
    annual_starter_rate = zeros(11, prim.n_sims)
    annual_stopper_rate = zeros(11, prim.n_sims)
    # Add 1 and 0 to export decisions once to make sure neither are ever 0
    annual_firms_export_decisions_no_zero = hcat(annual_firms_export_decisions, ones(12, 1, prim.n_sims), zeros(12, 1, prim.n_sims))
    # Add 1 and -1 to export changes once to make sure neither are ever 0
    annual_firms_export_change_no_zero = hcat(annual_firms_export_change, ones(12, 1, prim.n_sims), -1 .* ones(12, 1, prim.n_sims)) 
    exporter_i = countmap(annual_firms_export_change_no_zero[1,:,:])
    exporter_i_minus1 = countmap(annual_firms_export_change_no_zero[1,:,:])

    for i = 2:12
        for k = 1:prim.n_sims
            exporter_i_minus1 = countmap(annual_firms_export_decisions_no_zero[i-1,:,k])
            exporter_i = countmap(annual_firms_export_change_no_zero[i,:,k])
            
            annual_starter_rate[i-1,k] = exporter_i[1]/exporter_i_minus1[0]
            annual_stopper_rate[i-1,k] = exporter_i[-1]/exporter_i_minus1[1]
        end
    end

    #########################################################################
    # Find the coefficient of variation of log domestic sales for every sim #
    #########################################################################
    coef_variation = zeros(prim.n_sims)
    for i = 1:prim.n_sims
        coef_variation[i] = std(log.(annual_firms_sales_domestic[:,:,i]))/mean(log.(annual_firms_sales_domestic[:,:,i]))
    end

    ##############################################
    # Find the coefficient on log domestic sales #
    ##############################################
    β_sales = zeros(1)
    for s = 1:1
        # sales_today = mean(annual_firms_sales[2,:,:], dims=2)
        # sales_yesterday = mean(annual_firms_sales[1,:,:], dims=2)
        sales_today = annual_firms_sales_domestic[2,:,s]
        sales_yesterday = annual_firms_sales_domestic[1,:,s]
        year = []
        firm_num = []
        for i = 3:12
            # sales_today = [sales_today; mean(annual_firms_sales[i,:,:],dims=2)]
            # sales_yesterday = [sales_yesterday; mean(annual_firms_sales[i-1,:,:],dims=2)]
            sales_today = [sales_today; annual_firms_sales_domestic[i,:,s]]
            sales_yesterday = [sales_yesterday; annual_firms_sales_domestic[i-1,:,s]]
        end

        for i = 2:12
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
        β_sales[s] = β_moment[1]
    end

    for k = 1:prim.n_sims
        n_no_entry = 0
        n_one_entry = 0
        n_two_entries = 0
        n_three_entries = 0
        n_four_entries = 0
        n_five_entries = 0
        n_six_entries = 0
        n_seven_entries = 0
        export_entries = Vector{Int64}()
        export_entries_no_zero = Vector{Int64}()

        for i = 1:prim.n_firms
            exporter_i_map = countmap(vcat(annual_firms_export_change[:,i,k],1))
            if exporter_i_map[1] == 2
                n_one_entry = n_one_entry + 1
                export_entries = push!(export_entries, 1)
                export_entries_no_zero = push!(export_entries_no_zero, 1)
            elseif exporter_i_map[1] == 3
                n_two_entries = n_two_entries + 1
                export_entries =  push!(export_entries, 2)
                export_entries_no_zero = push!(export_entries_no_zero, 2)
            elseif exporter_i_map[1] == 4
                n_three_entries = n_three_entries + 1
                export_entries =  push!(export_entries, 3)
                export_entries_no_zero = push!(export_entries_no_zero, 3)
            elseif exporter_i_map[1] == 5
                n_four_entries = n_four_entries + 1
                export_entries =  push!(export_entries, 4)
                export_entries_no_zero = push!(export_entries_no_zero, 4)
            elseif exporter_i_map[1] == 6
                n_five_entries = n_five_entries + 1
                export_entries =  push!(export_entries, 5)
                export_entries_no_zero = push!(export_entries_no_zero, 5)
            elseif exporter_i_map[1] == 7
                n_six_entries = n_six_entries + 1
                export_entries =  push!(export_entries, 6)
                export_entries_no_zero = push!(export_entries_no_zero, 6)
            elseif exporter_i_map[1] == 8
                n_seven_entries = n_seven_entries + 1
                export_entries =  push!(export_entries, 7)
                export_entries_no_zero = push!(export_entries_no_zero, 7)
            elseif exporter_i_map[1] == 1
                n_no_entry = n_no_entry + 1
                export_entries =  push!(export_entries, 0)
            end
        end

        prop_two_entry[k] = n_two_entries/(prim.n_firms-n_no_entry)
        prop_three_entry[k] = n_three_entries/(prim.n_firms-n_no_entry)
        prop_four_entry[k] = n_four_entries/(prim.n_firms-n_no_entry)
        prop_five_entry[k] = n_five_entries/(prim.n_firms-n_no_entry)
        prop_six_entry[k] = n_six_entries/(prim.n_firms-n_no_entry)
        prop_seven_entry[k] = n_seven_entries/(prim.n_firms-n_no_entry)
        prop_one_entry[k] = n_one_entry/(prim.n_firms-n_no_entry)
        
        prop_no_entry_all[k] = n_no_entry/(prim.n_firms)
        prop_one_entry_all[k] = n_one_entry/(prim.n_firms)
        prop_two_entry_all[k] = n_two_entries/(prim.n_firms)
        prop_three_entry_all[k] = n_three_entries/(prim.n_firms)
        prop_four_entry_all[k] = n_four_entries/(prim.n_firms)
    end

    prop_one_entry_mean = round(mean(prop_one_entry), digits = 4)
    prop_two_entry_mean = round(mean(prop_two_entry), digits = 4)
    prop_three_entry_mean = round(mean(prop_three_entry), digits = 4)
    prop_four_entry_mean = round(mean(prop_four_entry), digits = 4)
    prop_five_entry_mean = round(mean(prop_five_entry), digits = 4)
    
    prop_no_entry_all_mean = round(mean(prop_no_entry_all), digits = 4)
    prop_one_entry_all_mean = round(mean(prop_one_entry_all), digits = 4)
    prop_two_entry_all_mean = round(mean(prop_two_entry_all), digits = 4)
    prop_three_entry_all_mean = round(mean(prop_three_entry_all), digits = 4)
    prop_four_entry_all_mean = round(mean(prop_four_entry_all), digits = 4)

    # distribution_SSE = (prop_two_entry_mean-prim.true_enter_twice_perc)^2+(prop_one_entry_mean-prim.true_enter_once_perc)^2+(prop_three_entry_mean-prim.true_enter_thrice_perc)^2+(prop_four_entry_mean-prim.true_enter_four_perc)^2+(prop_five_entry_mean-prim.true_enter_five_perc)^2

    exits_all = Int.(annual_firms_export_change .== -1)
    exits_cumulative = Int.(annual_firms_export_change .== -1)
    final_exit_no_reentry = ones(Int, 12, prim.n_firms, prim.n_sims)
    entries_all = Int.(annual_firms_export_change .== 1)

    for i in 2:12
        exits_cumulative[i,:,:] .+= exits_cumulative[i-1,:,:]
    end

    ##############################################
    # Calculate average years out before reentry #
    ##############################################
    # Set entries where the firm is exporting equal to zero, keeping strings of numbers with length equal to years out of each exit
    exits_cumulative = exits_cumulative .* (ones(Int, 12, prim.n_firms, prim.n_sims) .- annual_firms_export_decisions[:,:,:])
    # Test if that exit is the final exit and has no re-entry tied to it
    for i in 1:11
        final_exit_no_reentry[i,:,:] = Int.((exits_cumulative[i,:,:] .== exits_cumulative[12,:,:]))
    end
    # Only keep exits that result in re-entries
    exits_cumulative = exits_cumulative .* (1 .- final_exit_no_reentry)
    # Loop through the number of exits and sum up the years out that result in re-entry
    max_exits = maximum(exits_cumulative)
    exits_count_map = countmap(exits_cumulative)
    years_out_before_reentry = 0
    for i in 1:max_exits
        # Test if that number of exits actually occurs ever
        try
            exits_count_map[i]
        catch
            # if not set it to 0
            exits_count_map[i] = 0
        end
        years_out_before_reentry += exits_count_map[i]
    end
    # Get maximum exits for each firm-simulation, which is number of re-entries
    reentries = sum(maximum(exits_cumulative, dims = 1))
    # If there are no re-entries make them very small so it will be very far off the moment
    if reentries == 0
        reentries = 1e-5
        years_out_before_reentry = 1
    end

    ################################################################
    # Calculate percentage of exiters that re-enter after one year #
    ################################################################
    exits_all_noEndYears = sum(exits_all[2:11, :, :])
    immediate_reentries = sum(Int.((exits_all[2:11, :, :] .== 1) .& (entries_all[3:12, :, :] .== 1)))

    # If there are no exits, change it such that the data will be very far off the moment
    if exits_all_noEndYears == 0
        exits_all_noEndYears = 1e-5
        immediate_reentries = 1
    end

    ###########################
    # Calculate output vector #
    ###########################
    output = zeros(11)

    output[1] = mean(annual_starter_rate)
    output[2] = mean(annual_stopper_rate)
    output[3] = sum(annual_firms_sales_foreign./annual_firms_sales)/sum(annual_firms_export_decisions)
    output[4] = mean(coef_variation)
    # output[5] = β_moment[1]
    output[5] = mean(β_sales)
    output[6] = years_out_before_reentry/reentries
    output[7] = immediate_reentries/exits_all_noEndYears
    output[8] = prop_one_entry_mean
    output[9] = prop_two_entry_mean
    output[10] = prop_three_entry_mean
    output[11] = prop_four_entry_mean

    if model == 1
        error = abs(output[1]-prim.true_starter)/prim.true_starter+abs(output[2]-prim.true_stopper)/prim.true_stopper + abs(output[3]-prim.true_ave_es_ratio)/prim.true_ave_es_ratio+abs.(output[4]-prim.true_coef_var)/prim.true_coef_var+abs.(output[5]-prim.true_a_exp_growth)/prim.true_a_exp_growth 
    else
        error = abs(output[1]-prim.true_starter)/prim.true_starter+abs(output[2]-prim.true_stopper)/prim.true_stopper+abs(output[6]-prim.true_avg_time_out_reentrant)/prim.true_avg_time_out_reentrant+abs(output[7]-prim.true_perc_reenter_immediate)/prim.true_perc_reenter_immediate + abs(output[3]-prim.true_ave_es_ratio)/prim.true_ave_es_ratio+abs.(output[4]-prim.true_coef_var)/prim.true_coef_var+abs.(output[5]-prim.true_a_exp_growth)/prim.true_a_exp_growth 
    end
    
    if res.FC_0 < 0 || res.FC_1 < 0 || res.α_d < 0 || res.β_d < 0 || res.β_sq_d < 0 || res.δ < 0
        error = 100
    end
    
    print(error)
    print("\n")
    print(x)
    print("\n")
    print(output)
    print("\n")
    print([prim.true_starter, prim.true_stopper, prim.true_ave_es_ratio, prim.true_coef_var, prim.true_a_exp_growth, prim.true_avg_time_out_reentrant, prim.true_perc_reenter_immediate, prim.true_enter_once_perc, prim.true_enter_twice_perc, prim.true_enter_thrice_perc, prim.true_enter_four_perc]) 
    print("\n")
    print([prop_no_entry_all_mean, prop_one_entry_all_mean, prop_two_entry_all_mean, prop_three_entry_all_mean, prop_four_entry_all_mean])
    print("\n")
    print(median(annual_firms_sales))
    print("\n")
    print(median(annual_firms_sales_foreign))
    print("\n")
    print(mean(annual_firms_sales_foreign))
    print("\n")
    print(median(firms_labor_decisions[101:112,:,:]))
    print("\n")
    print(median(firms_capital_decisions[101:112,:,:]))
    print("\n")

    return error 

end

function data_sim_delta(prim::Primitives, res::Results)
    @unpack val_func, ex_func, n_func, k_func, n_prev_ex, ex_cap, ϵ, Q, σ_e, ρ_e, C_star = res #unpack value function
    @unpack Q_grid, ϵ_grid, n_periods, n_firms, ρ_q, σ_q, θ, α_n = prim #unpack primitives

    Random.seed!(17)

    firms_export_capital = ones(n_periods, n_firms)
    firms_export_decisions = zeros(n_periods, n_firms)
    firms_labor_decisions = ones(n_periods, n_firms)
    firms_capital_decisions = ones(n_periods, n_firms)
    firms_sales_domestic = zeros(n_periods, n_firms)
    firms_sales = zeros(n_periods, n_firms)

    for i = 2:n_periods

        Q[i] = exp(ρ_q*log(Q[i-1]) + rand(Normal(0, σ_q)))
        if Q[i] > 1.2
            Q[i] == 1.2
        elseif Q[i] < 0.8
            Q[i] == 0.8
        end
        Q_index = findmin(abs.(Q[i] .- Q_grid))[2]

        for j = 1:n_firms
            
            ϵ[i,j] = exp(ρ_e*log(ϵ[i-1,j]) + rand(Normal(0,σ_e)))
            ϵ_index = findmin(abs.(ϵ[i,j] .- ϵ_grid))[2]

            firms_export_decisions[i,j] = ex_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j])]
            if firms_export_decisions[i,j] == 0 && firms_export_capital[i-1,j] > 1
                firms_export_capital[i,j] = firms_export_capital[i-1,j] - 1
            elseif firms_export_decisions[i,j] == 1
                firms_export_capital[i,j] = n_prev_ex
            end
            firms_labor_decisions[i,j] = n_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j])]
            firms_capital_decisions[i,j] = k_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j])]
            firms_sales_domestic[i,j] = val_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j])] + prim.w*n_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j])] + prim.r*k_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j])]-firms_export_decisions[i,j]*(Q_grid[Q_index]*(1+res.τ))^θ*C_star^(1/θ)*ϵ_grid[ϵ_index]*firms_labor_decisions[i,j]^(α_n*(θ-1)/θ)*firms_capital_decisions[i,j]^((1-α_n)*(θ-1)/θ)
            firms_sales[i,j] = val_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i,j])] + prim.w*n_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i,j])] + prim.r*k_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i,j])]
        end
    end

    return firms_export_decisions, firms_labor_decisions, firms_capital_decisions, firms_sales_domestic, firms_sales
end

function profit_func(prim::Primitives, res::Results, x)
    # x[1] is labor choice, x[2] is capital, x[3] is export decision, x[4] is productivity, and x[5] is export capital, x[6] is Q, x[7] is FC_0, x[8] is FC_1
    @unpack θ, C, r, w, α_n = prim
    @unpack τ, C_star, FC_0, FC_1 = res
    Π = ((1+x[3]*(x[6]*(1-res.τ))^θ*C_star/C)^(1/θ)*C^(1/θ)*x[4]*x[1]^(α_n*(θ-1)/θ)*x[2]^((1-α_n)*(θ-1)/θ)-w*x[1]-r*x[2])-1*x[3]*((1-x[5])*x[7]+x[8])
    return Π
end

function export_revenue(prim::Primitives, res::Results, x)
    # x[1] is labor choice, x[2] is capital, x[3] is export decision, x[4] is productivity, and x[5] is Q
    @unpack θ, C, r, w, α_n = prim
    @unpack τ, C_star, FC_0, FC_1 = res
    # export_revenue = ((1+x[3]*(x[5]*(1-res.τ))^θ*C_star/C)^(1/θ)-1)*C^(1/θ)*x[4]*x[1]^(α_n*(θ-1)/θ)*x[2]^((1-α_n)*(θ-1)/θ)
    export_revenue = x[3]*x[5]*C_star^(1/θ)*(1/(1+x[5]^(-1*θ)*C/C_star)*x[4]^(θ/(θ-1))*x[1]^α_n*x[2]^(1-α_n))^((θ-1)/θ)
    return export_revenue
end

function domestic_revenue(prim::Primitives, res::Results, x)
    # x[1] is labor choice, x[2] is capital, x[3] is export decision, x[4] is productivity, and x[5] is Q
    @unpack θ, C, r, w, α_n = prim
    @unpack τ, C_star, FC_0, FC_1 = res
    # domestic_revenue = ((1)^(1/θ)-1)*C^(1/θ)*x[4]*x[1]^(α_n*(θ-1)/θ)*x[2]^((1-α_n)*(θ-1)/θ)
    domestic_revenue = C^(1/θ)*((1-x[3]*1/(1+x[5]^(-1*θ)*C/C_star))*x[4]^(θ/(θ-1))*x[1]^α_n*x[2]^(1-α_n))^((θ-1)/θ)
    return domestic_revenue
end

function total_revenue(prim::Primitives, res::Results, x)
    # x[1] is labor choice, x[2] is capital, x[3] is export decision, x[4] is productivity, and x[5] is Q
    @unpack θ, C, r, w, α_n = prim
    @unpack τ, C_star, FC_0, FC_1 = res
    # total_revenue = (1+x[3]*(x[5]*(1-res.τ))^θ*C_star/C)^(1/θ)*C^(1/θ)*x[4]*x[1]^(α_n*(θ-1)/θ)*x[2]^((1-α_n)*(θ-1)/θ)
    total_revenue = domestic_revenue(prim, res, x) + export_revenue(prim, res, x)
    return total_revenue
end

function tariff_experiment(prim::Primitives, res::Results, tariff::Int64, solve::Int64, filename::AbstractString)
    Random.seed!(17)
    @unpack Q_grid, ϵ_grid, n_periods_experiment, n_firms, n_sims, ρ_q, σ_q = prim #unpack primitives
    @unpack Q_experiment, ϵ_experiment, n_prev_ex, ex_cap, ρ_e, σ_e = res #unpack results

    if solve == 1
        res.τ = 0/100
        Solve_model(prim, res)
        
        save_object("./model/objects/normal_val_func_$filename.jld2", res.val_func)
        save_object("./model/objects/normal_ex_func_$filename.jld2", res.ex_func)
        save_object("./model/objects/normal_n_func_$filename.jld2", res.n_func)
        save_object("./model/objects/normal_k_func_$filename.jld2", res.k_func)

        normal_val_func = copy(res.val_func)
        normal_ex_func = copy(res.ex_func)
        normal_n_func = copy(res.n_func)
        normal_k_func = copy(res.k_func)
        
        if tariff > 0
            res.τ = tariff/100
            Solve_model(prim, res)
            tariff_val_func = copy(res.val_func)
            tariff_ex_func = copy(res.ex_func)
            tariff_n_func = copy(res.n_func)
            tariff_k_func = copy(res.k_func)
            
            save_object("./model/objects/tariff_val_func_$filename.jld2", res.val_func)
            save_object("./model/objects/tariff_ex_func_$filename.jld2", res.ex_func)
            save_object("./model/objects/tariff_n_func_$filename.jld2", res.n_func)
            save_object("./model/objects/tariff_k_func_$filename.jld2", res.k_func)
    
        else
            res.τ = 0/100
            tariff_val_func = copy(res.val_func)
            tariff_ex_func = copy(res.ex_func)
            tariff_n_func = copy(res.n_func)
            tariff_k_func = copy(res.k_func)
        end
    else
        normal_val_func = load_object("./model/objects/normal_val_func_$filename.jld2")
        normal_ex_func = load_object("./model/objects/normal_ex_func_$filename.jld2")
        normal_n_func = load_object("./model/objects/normal_n_func_$filename.jld2")
        normal_k_func = load_object("./model/objects/normal_k_func_$filename.jld2")
        
        tariff_val_func = load_object("./model/objects/tariff_val_func_$filename.jld2")
        tariff_ex_func = load_object("./model/objects/tariff_ex_func_$filename.jld2")
        tariff_n_func = load_object("./model/objects/tariff_n_func_$filename.jld2")
        tariff_k_func = load_object("./model/objects/tariff_k_func_$filename.jld2")
    end

    firms_export_capital = ones(n_periods_experiment, n_firms, n_sims)
    firms_export_decisions = zeros(n_periods_experiment, n_firms, n_sims)
    firms_labor_decisions = ones(n_periods_experiment, n_firms, n_sims)
    firms_capital_decisions = ones(n_periods_experiment, n_firms, n_sims)
    firms_sales_non_exporter = zeros(n_periods_experiment, n_firms, n_sims)
    firms_sales = zeros(n_periods_experiment, n_firms, n_sims)
    firms_export_sales = zeros(n_periods_experiment, n_firms, n_sims)

    for k = 1:n_sims
        Random.seed!(k)
        for i = 2:n_periods_experiment

            Q_experiment[i] = exp(ρ_q*log(Q_experiment[i-1]) + rand(Normal(0, σ_q)))
            if Q_experiment[i] > 3
                Q_experiment[i] == 3
            end
            Q_index = findmin(abs.(Q_experiment[i] .- Q_grid))[2]
            Q_experiment[i] = Q_grid[Q_index]

            for j = 1:n_firms
                
                ϵ_experiment[i,j] = exp(ρ_e*log(ϵ_experiment[i-1,j]) + rand(Normal(0,σ_e)))
                ϵ_index = findmin(abs.(ϵ_experiment[i,j] .- ϵ_grid))[2]
                ϵ_experiment[i,j] = ϵ_grid[ϵ_index]

                if i == 106
                    firms_export_decisions[i,j,k] = tariff_ex_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                else
                    firms_export_decisions[i,j,k] = normal_ex_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                end
                
                if firms_export_decisions[i,j,k] == 0 && firms_export_capital[i-1,j,k] > 1
                    firms_export_capital[i,j,k] = firms_export_capital[i-1,j,k] - 1
                elseif firms_export_decisions[i,j,k] == 1
                    firms_export_capital[i,j,k] = n_prev_ex
                end
                
                if i == 106
                    firms_labor_decisions[i,j,k] = tariff_n_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                    firms_capital_decisions[i,j,k] = tariff_k_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                    res.τ = tariff/100
                else
                    firms_labor_decisions[i,j,k] = normal_n_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                    firms_capital_decisions[i,j,k] = normal_k_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                    res.τ = 0
                end
                firms_sales_non_exporter[i,j,k] = domestic_revenue(prim, res, [firms_labor_decisions[i,j,k] firms_capital_decisions[i,j,k] firms_export_decisions[i,j,k] ϵ_grid[ϵ_index] Q_grid[Q_index]])
                firms_sales[i,j,k] = total_revenue(prim, res, [firms_labor_decisions[i,j,k] firms_capital_decisions[i,j,k] firms_export_decisions[i,j,k] ϵ_grid[ϵ_index] Q_grid[Q_index]])
                firms_export_sales[i,j,k] = export_revenue(prim, res, [firms_labor_decisions[i,j,k] firms_capital_decisions[i,j,k] firms_export_decisions[i,j,k] ϵ_grid[ϵ_index] Q_grid[Q_index]])
            end
        end
    end

    return firms_export_decisions, firms_labor_decisions, firms_capital_decisions, firms_sales_non_exporter, firms_sales, firms_export_sales

end

function Q_experiment(prim::Primitives, res::Results, filename::AbstractString)
    Random.seed!(17)
    @unpack Q_grid, ϵ_grid, n_periods_experiment, n_firms, n_sims, ρ_q, σ_q = prim #unpack primitives
    @unpack Q_experiment, ϵ_experiment, n_prev_ex, ex_cap, ρ_e, σ_e = res #unpack results

    normal_val_func = load_object("./model/objects/normal_val_func_$filename.jld2")
    normal_ex_func = load_object("./model/objects/normal_ex_func_$filename.jld2")
    normal_n_func = load_object("./model/objects/normal_n_func_$filename.jld2")
    normal_k_func = load_object("./model/objects/normal_k_func_$filename.jld2")

    firms_export_capital = ones(n_periods_experiment, n_firms, n_sims)
    firms_export_decisions = zeros(n_periods_experiment, n_firms, n_sims)
    firms_labor_decisions = ones(n_periods_experiment, n_firms, n_sims)
    firms_capital_decisions = ones(n_periods_experiment, n_firms, n_sims)
    firms_sales_non_exporter = zeros(n_periods_experiment, n_firms, n_sims)
    firms_sales = zeros(n_periods_experiment, n_firms, n_sims)

    for k = 1:n_sims
        Random.seed!(k)
        for i = 2:n_periods_experiment

            Q_index = 3
            
            if i < 105 || i > 107
                Q_experiment[i] = exp(ρ_q*log(Q_experiment[i-1]) + rand(Normal(0, σ_q)))
                if Q_experiment[i] > 3
                    Q_experiment[i] == 3
                end
                Q_index = findmin(abs.(Q_experiment[i] .- Q_grid))[2]
                Q_experiment[i] = Q_grid[Q_index]
            elseif i == 105 || i == 107
                Q_experiment[i] = 1
                Q_index = 3
            elseif i == 106
                Q_experiment[i] = 0.9
                Q_index = 2
            end

            for j = 1:n_firms
                
                ϵ_experiment[i,j] = exp(ρ_e*log(ϵ_experiment[i-1,j]) + rand(Normal(0,σ_e)))
                ϵ_index = findmin(abs.(ϵ_experiment[i,j] .- ϵ_grid))[2]
                
                firms_export_decisions[i,j,k] = normal_ex_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                if firms_export_decisions[i,j,k] == 0 && firms_export_capital[i-1,j,k] > 1
                    firms_export_capital[i,j,k] = firms_export_capital[i-1,j,k] - 1
                elseif firms_export_decisions[i,j,k] == 1
                    firms_export_capital[i,j,k] = n_prev_ex
                end
                
                firms_labor_decisions[i,j,k] = normal_n_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                firms_sales_non_exporter[i,j,k] = normal_val_func[Q_index, ϵ_index, 1] + prim.w*normal_n_func[Q_index, ϵ_index, 1] + prim.r*normal_k_func[Q_index, ϵ_index, 1]
                firms_sales[i,j,k] = normal_val_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i,j,k])] + prim.w*normal_n_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i,j,k])] + prim.r*normal_k_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i,j,k])]
            end
        end
    end

    return firms_export_decisions, firms_labor_decisions, firms_capital_decisions, firms_sales_non_exporter, firms_sales

end

function export_experience_experiment(prim::Primitives, res::Results, filename::AbstractString)
    Random.seed!(17)
    @unpack Q_grid, ϵ_grid, n_periods_experiment, n_firms, n_sims, ρ_q, σ_q = prim #unpack primitives
    @unpack Q_experiment, ϵ_experiment, n_prev_ex, ex_cap, ρ_e, σ_e = res #unpack results

    normal_val_func = load_object("./model/objects/normal_val_func_$filename.jld2")
    normal_ex_func = load_object("./model/objects/normal_ex_func_$filename.jld2")
    normal_n_func = load_object("./model/objects/normal_n_func_$filename.jld2")
    normal_k_func = load_object("./model/objects/normal_k_func_$filename.jld2")

    firms_export_capital = ones(n_periods_experiment, n_firms, n_sims)
    firms_export_decisions = zeros(n_periods_experiment, n_firms, n_sims)
    firms_labor_decisions = ones(n_periods_experiment, n_firms, n_sims)
    firms_capital_decisions = ones(n_periods_experiment, n_firms, n_sims)
    firms_sales_non_exporter = zeros(n_periods_experiment, n_firms, n_sims)
    firms_sales = zeros(n_periods_experiment, n_firms, n_sims)
    firms_export_sales = zeros(n_periods_experiment, n_firms, n_sims)

    for k = 1:n_sims
        Random.seed!(k)
        for i = 2:n_periods_experiment
            
            Q_experiment[i] = exp(ρ_q*log(Q_experiment[i-1]) + rand(Normal(0, σ_q)))
            if Q_experiment[i] > 3
                Q_experiment[i] == 3
            end
            Q_index = findmin(abs.(Q_experiment[i] .- Q_grid))[2]
            Q_experiment[i] = Q_grid[Q_index]

            for j = 1:n_firms
                
                ϵ_experiment[i,j] = exp(ρ_e*log(ϵ_experiment[i-1,j]) + rand(Normal(0,σ_e)))
                ϵ_index = findmin(abs.(ϵ_experiment[i,j] .- ϵ_grid))[2]
                ϵ_experiment[i,j] = ϵ_grid[ϵ_index]
                
                firms_export_decisions[i,j,k] = normal_ex_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                if firms_export_decisions[i,j,k] == 0 && firms_export_capital[i-1,j,k] > 1
                    firms_export_capital[i,j,k] = firms_export_capital[i-1,j,k] - 1
                elseif firms_export_decisions[i,j,k] == 1
                    firms_export_capital[i,j,k] = n_prev_ex
                end
                
                firms_labor_decisions[i,j,k] = normal_n_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                firms_sales_non_exporter[i,j,k] = normal_val_func[Q_index, ϵ_index, 1] + prim.w*normal_n_func[Q_index, ϵ_index, 1] + prim.r*normal_k_func[Q_index, ϵ_index, 1]
                firms_sales[i,j,k] = normal_val_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i,j,k])] + prim.w*normal_n_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i,j,k])] + prim.r*normal_k_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i,j,k])]
                firms_export_sales[i,j,k] = export_revenue(prim, res, [firms_labor_decisions[i,j,k] firms_capital_decisions[i,j,k] firms_export_decisions[i,j,k] ϵ_experiment[i,j] Q_experiment[i]])

                # After choices are made, export capital of everyone decays to 0 in period 106 and we continue
                if i == 106
                    firms_export_capital[i,j,k] = 1
                end
            end
        end
    end

    return firms_export_decisions, firms_labor_decisions, firms_capital_decisions, firms_sales_non_exporter, firms_sales, firms_export_sales

end

function export_subsidy_experiment(prim::Primitives, res::Results)
    Random.seed!(17)
    # But choice to export does depend on continuation value... Makes one time shock hard
    # Use continuation value from the standard model solution but the profits today be with FC_0 = 0 ??

    # Initialize to set standard parameters
    prim, res = Initialize(3)
    subsidy_rate = 0.5

    Solve_model(prim, res)
    # Doing subsidy afterwards makes sure continuation value is accurate
    V_iterate_subsidy(prim, res, subsidy_rate)

    @unpack Q_grid, ϵ_grid, n_periods_experiment, n_firms, n_sims, ρ_q, σ_q = prim #unpack primitives
    @unpack Q_experiment, ϵ_experiment, n_prev_ex, prev_ex_grid, ex_cap, ρ_e, σ_e, FC_0, FC_1 = res #unpack results

    firms_export_capital_rand = ones(n_periods_experiment, n_firms, n_sims)
    firms_export_decisions_rand = zeros(n_periods_experiment, n_firms, n_sims)
    firms_labor_decisions_rand = ones(n_periods_experiment, n_firms, n_sims)
    firms_capital_decisions_rand = ones(n_periods_experiment, n_firms, n_sims)
    firms_sales_non_exporter_rand = zeros(n_periods_experiment, n_firms, n_sims)
    firms_sales_rand = zeros(n_periods_experiment, n_firms, n_sims)
    firms_export_sales_rand = zeros(n_periods_experiment, n_firms, n_sims)
    
    firms_export_capital_targeted = ones(n_periods_experiment, n_firms, n_sims)
    firms_export_decisions_targeted = zeros(n_periods_experiment, n_firms, n_sims)
    firms_labor_decisions_targeted = ones(n_periods_experiment, n_firms, n_sims)
    firms_capital_decisions_targeted = ones(n_periods_experiment, n_firms, n_sims)
    firms_sales_non_exporter_targeted = zeros(n_periods_experiment, n_firms, n_sims)
    firms_sales_targeted = zeros(n_periods_experiment, n_firms, n_sims)
    firms_export_sales_targeted = zeros(n_periods_experiment, n_firms, n_sims)

    productivity_exports = zeros(n_firms, 3, n_sims)

    rand_firm = shuffle(Xoshiro(17), collect(1:n_firms))

    for k = 1:n_sims
        Random.seed!(k)
        subsidy_budget_rand = 0
        n_subsidy_rand = 0
        subsidy_budget_targeted = 0
        n_subsidy_targeted = 0
        
        for i = 2:n_periods_experiment
            
            if i == 105 || i == 106
                Q_experiment[i] = 0.8 #exp(ρ_q*log(Q_experiment[i-1]) + rand(Normal(0, σ_q))) # 1
            else
                Q_experiment[i] = 0.8 #exp(ρ_q*log(Q_experiment[i-1]) + rand(Normal(0, σ_q)))
            end
            Q_index = findmin(abs.(Q_experiment[i] .- Q_grid))[2]

            # Update ϵ's first:
            for j = 1:n_firms
                # Update ϵ for non-subsidy year
                ϵ_experiment[i,j] = exp(ρ_e*log(ϵ_experiment[i-1,j]) + rand(Normal(0,σ_e)))
                ϵ_index = findmin(abs.(ϵ_experiment[i,j] .- ϵ_grid))[2]
                ϵ_experiment[i,j] = ϵ_grid[ϵ_index]
            end

            if i == 106
                target_firm = sortperm(firms_export_capital_targeted[i-1,:,k], rev=true)
                productivity_exports[:,1,k] = transpose(ϵ_experiment[i,:])
            end

            for j = 1:n_firms

                if i != 106
                    # Normal times:
                    ϵ_index = findmin(abs.(ϵ_experiment[i,j] .- ϵ_grid))[2]
                    firms_export_decisions_rand[i,j,k] = res.ex_func[Q_index, ϵ_index, floor(Int,firms_export_capital_rand[i-1,j,k])]
                    if firms_export_decisions_rand[i,j,k] == 0 && firms_export_capital_rand[i-1,j,k] > 1
                        firms_export_capital_rand[i,j,k] = firms_export_capital_rand[i-1,j,k] - 1
                    elseif firms_export_decisions_rand[i,j,k] == 1
                        firms_export_capital_rand[i,j,k] = n_prev_ex
                    end
                    
                    firms_labor_decisions_rand[i,j,k] = res.n_func[Q_index, ϵ_index, floor(Int,firms_export_capital_rand[i-1,j,k])]
                    firms_capital_decisions_rand[i,j,k] = res.k_func[Q_index, ϵ_index, floor(Int,firms_export_capital_rand[i-1,j,k])]
                    firms_sales_non_exporter_rand[i,j,k] = domestic_revenue(prim, res, [firms_labor_decisions_rand[i,j,k] firms_capital_decisions_rand[i,j,k] firms_export_decisions_rand[i,j,k] ϵ_experiment[i,j] Q_experiment[i]])
                    firms_sales_rand[i,j,k] = total_revenue(prim, res, [firms_labor_decisions_rand[i,j,k] firms_capital_decisions_rand[i,j,k] firms_export_decisions_rand[i,j,k] ϵ_experiment[i,j] Q_experiment[i]])
                    firms_export_sales_rand[i,j,k] = export_revenue(prim, res, [firms_labor_decisions_rand[i,j,k] firms_capital_decisions_rand[i,j,k] firms_export_decisions_rand[i,j,k] ϵ_experiment[i,j] Q_experiment[i]])

                    firms_export_decisions_targeted[i,j,k] = res.ex_func[Q_index, ϵ_index, floor(Int,firms_export_capital_targeted[i-1,j,k])]
                    if firms_export_decisions_targeted[i,j,k] == 0 && firms_export_capital_targeted[i-1,j,k] > 1
                        firms_export_capital_targeted[i,j,k] = firms_export_capital_targeted[i-1,j,k] - 1
                    elseif firms_export_decisions_targeted[i,j,k] == 1
                        firms_export_capital_targeted[i,j,k] = n_prev_ex
                    end
                    
                    firms_labor_decisions_targeted[i,j,k] = res.n_func[Q_index, ϵ_index, floor(Int,firms_export_capital_targeted[i-1,j,k])]
                    firms_capital_decisions_targeted[i,j,k] = res.k_func[Q_index, ϵ_index, floor(Int,firms_export_capital_targeted[i-1,j,k])]
                    firms_sales_non_exporter_targeted[i,j,k] = domestic_revenue(prim, res, [firms_labor_decisions_targeted[i,j,k] firms_capital_decisions_targeted[i,j,k] firms_export_decisions_targeted[i,j,k] ϵ_experiment[i,j] Q_experiment[i]])
                    firms_sales_targeted[i,j,k] = total_revenue(prim, res, [firms_labor_decisions_targeted[i,j,k] firms_capital_decisions_targeted[i,j,k] firms_export_decisions_targeted[i,j,k] ϵ_experiment[i,j] Q_experiment[i]])
                    firms_export_sales_targeted[i,j,k] = export_revenue(prim, res, [firms_labor_decisions_targeted[i,j,k] firms_capital_decisions_targeted[i,j,k] firms_export_decisions_targeted[i,j,k] ϵ_experiment[i,j] Q_experiment[i]])
                    if i == 105
                        subsidy_budget_rand += profit_func(prim, res, [firms_labor_decisions_rand[i,j,k], firms_capital_decisions_rand[i,j,k], 
                                                                       firms_export_decisions_rand[i,j,k], ϵ_experiment[i,j], 
                                                                       firms_export_capital_rand[i,j,k], Q_experiment[i], FC_0, FC_1])*0.01
                        subsidy_budget_targeted += profit_func(prim, res, [firms_labor_decisions_targeted[i,j,k], firms_capital_decisions_targeted[i,j,k], 
                                                                       firms_export_decisions_targeted[i,j,k], ϵ_experiment[i,j],
                                                                       firms_export_capital_targeted[i,j,k], Q_experiment[i], FC_0, FC_1])*0.01
                    end            
                else
                    # UPDATE INDICES FOR RAND AND SUBSIDY
                    ϵ_index = findmin(abs.(ϵ_experiment[i,rand_firm[j]] .- ϵ_grid))[2]
                    if subsidy_budget_rand > 0
                        # Get's subsidy if any is available. If export capital = 1, subsidy = 0
                        firms_export_decisions_rand[i,rand_firm[j],k] = res.ex_func_subsidy[Q_index, ϵ_index, floor(Int,firms_export_capital_rand[i-1,rand_firm[j],k])]
                        if firms_export_decisions_rand[i,rand_firm[j],k] == 0 && firms_export_capital_rand[i-1,rand_firm[j],k] > 1
                            firms_export_capital_rand[i,rand_firm[j],k] = firms_export_capital_rand[i-1,rand_firm[j],k] - 1
                        elseif firms_export_decisions_rand[i,rand_firm[j],k] == 1
                            firms_export_capital_rand[i,rand_firm[j],k] = n_prev_ex
                        end
                        
                        firms_labor_decisions_rand[i,rand_firm[j],k] = res.n_func_subsidy[Q_index, ϵ_index, floor(Int,firms_export_capital_rand[i-1,rand_firm[j],k])]
                        firms_capital_decisions_rand[i,rand_firm[j],k] = res.k_func_subsidy[Q_index, ϵ_index, floor(Int,firms_export_capital_rand[i-1,rand_firm[j],k])]
                        firms_sales_non_exporter_rand[i,rand_firm[j],k] = domestic_revenue(prim, res, [firms_labor_decisions_rand[i,rand_firm[j],k] firms_capital_decisions_rand[i,rand_firm[j],k] firms_export_decisions_rand[i,rand_firm[j],k] ϵ_experiment[i,rand_firm[j]] Q_experiment[i]])
                        firms_sales_rand[i,rand_firm[j],k] = total_revenue(prim, res, [firms_labor_decisions_rand[i,rand_firm[j],k] firms_capital_decisions_rand[i,rand_firm[j],k] firms_export_decisions_rand[i,rand_firm[j],k] ϵ_experiment[i,rand_firm[j]] Q_experiment[i]])
                        firms_export_sales_rand[i,rand_firm[j],k] = export_revenue(prim, res, [firms_labor_decisions_rand[i,rand_firm[j],k] firms_capital_decisions_rand[i,rand_firm[j],k] firms_export_decisions_rand[i,rand_firm[j],k] ϵ_experiment[i,rand_firm[j]] Q_experiment[i]])

                        # Subtract the cost from the budget and add one to counter
                        n_subsidy_rand += 1
                        subsidy_budget_rand -= min((1-prev_ex_grid[floor(Int,firms_export_capital_rand[i-1,rand_firm[j],k])])*FC_0, subsidy_rate*FC_0)
                    else
                        firms_export_decisions_rand[i,rand_firm[j],k] = res.ex_func[Q_index, ϵ_index, floor(Int,firms_export_capital_rand[i-1,rand_firm[j],k])]
                        if firms_export_decisions_rand[i,rand_firm[j],k] == 0 && firms_export_capital_rand[i-1,rand_firm[j],k] > 1
                            firms_export_capital_rand[i,rand_firm[j],k] = firms_export_capital_rand[i-1,rand_firm[j],k] - 1
                        elseif firms_export_decisions_rand[i,rand_firm[j],k] == 1
                            firms_export_capital_rand[i,rand_firm[j],k] = n_prev_ex
                        end
                        
                        firms_labor_decisions_rand[i,rand_firm[j],k] = res.n_func[Q_index, ϵ_index, floor(Int,firms_export_capital_rand[i-1,rand_firm[j],k])]
                        firms_capital_decisions_rand[i,rand_firm[j],k] = res.k_func[Q_index, ϵ_index, floor(Int,firms_export_capital_rand[i-1,rand_firm[j],k])]
                        firms_sales_non_exporter_rand[i,rand_firm[j],k] = domestic_revenue(prim, res, [firms_labor_decisions_rand[i,rand_firm[j],k] firms_capital_decisions_rand[i,rand_firm[j],k] firms_export_decisions_rand[i,rand_firm[j],k] ϵ_experiment[i,rand_firm[j]] Q_experiment[i]])
                        firms_sales_rand[i,rand_firm[j],k] = total_revenue(prim, res, [firms_labor_decisions_rand[i,rand_firm[j],k] firms_capital_decisions_rand[i,rand_firm[j],k] firms_export_decisions_rand[i,rand_firm[j],k] ϵ_experiment[i,rand_firm[j]] Q_experiment[i]])
                        firms_export_sales_rand[i,rand_firm[j],k] = export_revenue(prim, res, [firms_labor_decisions_rand[i,rand_firm[j],k] firms_capital_decisions_rand[i,rand_firm[j],k] firms_export_decisions_rand[i,rand_firm[j],k] ϵ_experiment[i,rand_firm[j]] Q_experiment[i]])
                    end

                    ϵ_index = findmin(abs.(ϵ_experiment[i,target_firm[j]] .- ϵ_grid))[2]
                    if subsidy_budget_targeted > 0
                        # Get's subsidy if any is available
                        firms_export_decisions_targeted[i,target_firm[j],k] = res.ex_func_subsidy[Q_index, ϵ_index, floor(Int,firms_export_capital_targeted[i-1,target_firm[j],k])]
                        if firms_export_decisions_targeted[i,target_firm[j],k] == 0 && firms_export_capital_targeted[i-1,target_firm[j],k] > 1
                            firms_export_capital_targeted[i,target_firm[j],k] = firms_export_capital_targeted[i-1,target_firm[j],k] - 1
                        elseif firms_export_decisions_targeted[i,target_firm[j],k] == 1
                            firms_export_capital_targeted[i,target_firm[j],k] = n_prev_ex
                        end
                        
                        firms_labor_decisions_targeted[i,target_firm[j],k] = res.n_func_subsidy[Q_index, ϵ_index, floor(Int,firms_export_capital_targeted[i-1,target_firm[j],k])]
                        firms_capital_decisions_targeted[i,target_firm[j],k] = res.k_func_subsidy[Q_index, ϵ_index, floor(Int,firms_export_capital_rand[i-1,target_firm[j],k])]
                        firms_sales_non_exporter_targeted[i,target_firm[j],k] = domestic_revenue(prim, res, [firms_labor_decisions_targeted[i,target_firm[j],k] firms_capital_decisions_targeted[i,target_firm[j],k] firms_export_decisions_targeted[i,target_firm[j],k] ϵ_experiment[i,target_firm[j]] Q_experiment[i]])
                        firms_sales_targeted[i,target_firm[j],k] = total_revenue(prim, res, [firms_labor_decisions_targeted[i,target_firm[j],k] firms_capital_decisions_targeted[i,target_firm[j],k] firms_export_decisions_targeted[i,target_firm[j],k] ϵ_experiment[i,target_firm[j]] Q_experiment[i]])
                        firms_export_sales_targeted[i,target_firm[j],k] = export_revenue(prim, res, [firms_labor_decisions_targeted[i,target_firm[j],k] firms_capital_decisions_targeted[i,target_firm[j],k] firms_export_decisions_targeted[i,target_firm[j],k] ϵ_experiment[i,target_firm[j]] Q_experiment[i]])
                        # Subtract the cost from the budget and add one to counter
                        n_subsidy_targeted += 1
                        subsidy_budget_targeted -= min((1-prev_ex_grid[floor(Int,firms_export_capital_targeted[i-1,target_firm[j],k])])*FC_0, subsidy_rate*FC_0)
                    else
                        # Otherwise, no subsidy
                        firms_export_decisions_targeted[i,target_firm[j],k] = res.ex_func[Q_index, ϵ_index, floor(Int,firms_export_capital_targeted[i-1,target_firm[j],k])]
                        if firms_export_decisions_targeted[i,target_firm[j],k] == 0 && firms_export_capital_targeted[i-1,target_firm[j],k] > 1
                            firms_export_capital_targeted[i,target_firm[j],k] = firms_export_capital_targeted[i-1,target_firm[j],k] - 1
                        elseif firms_export_decisions_targeted[i,target_firm[j],k] == 1
                            firms_export_capital_targeted[i,target_firm[j],k] = n_prev_ex
                        end
                        
                        firms_labor_decisions_targeted[i,target_firm[j],k] = res.n_func[Q_index, ϵ_index, floor(Int,firms_export_capital_targeted[i-1,target_firm[j],k])]
                        firms_capital_decisions_targeted[i,target_firm[j],k] = res.k_func[Q_index, ϵ_index, floor(Int,firms_export_capital_targeted[i-1,target_firm[j],k])]
                        firms_sales_non_exporter_targeted[i,target_firm[j],k] = domestic_revenue(prim, res, [firms_labor_decisions_targeted[i,target_firm[j],k] firms_capital_decisions_targeted[i,target_firm[j],k] firms_export_decisions_targeted[i,target_firm[j],k] ϵ_experiment[i,target_firm[j]] Q_experiment[i]])
                        firms_sales_targeted[i,target_firm[j],k] = total_revenue(prim, res, [firms_labor_decisions_targeted[i,target_firm[j],k] firms_capital_decisions_targeted[i,target_firm[j],k] firms_export_decisions_targeted[i,target_firm[j],k] ϵ_experiment[i,target_firm[j]] Q_experiment[i]])
                        firms_export_sales_targeted[i,target_firm[j],k] = export_revenue(prim, res, [firms_labor_decisions_targeted[i,target_firm[j],k] firms_capital_decisions_targeted[i,target_firm[j],k] firms_export_decisions_targeted[i,target_firm[j],k] ϵ_experiment[i,target_firm[j]] Q_experiment[i]])
                    end
                    # somewhat wasteful with memory
                    productivity_exports[:,2,k] = transpose(firms_export_decisions_targeted[i,:,k])
                    productivity_exports[:,3,k] = transpose(firms_export_decisions_rand[i,:,k])
                end
            end
            # Print the budgets bc im curious!
            if i == 105 && k % 10 == 0
                print(subsidy_budget_rand)
                print("\n")
                print(subsidy_budget_targeted)
                print("\n")
            end
            # Print the budgets after bc im curious!
            if i == 107 && k % 10 == 0
                print(subsidy_budget_rand)
                print("\n")
                print(n_subsidy_rand)
                print("\n")
                print(subsidy_budget_targeted)
                print("\n")
                print(n_subsidy_targeted)
                print("\n")
            end
        end
    end

    return firms_export_decisions_rand, firms_labor_decisions_rand, firms_capital_decisions_rand, firms_sales_non_exporter_rand, firms_sales_rand, firms_export_sales_rand, firms_export_decisions_targeted, firms_labor_decisions_targeted, firms_capital_decisions_targeted, firms_sales_non_exporter_targeted, firms_sales_targeted, firms_export_sales_targeted, productivity_exports

end

function data_sim_delta_nsims(prim::Primitives, res::Results)
    @unpack val_func, ex_func, n_func, k_func, n_prev_ex, ex_cap, ϵ, Q, FC_0, FC_1, σ_e, ρ_e, C_star, tauchen_trans_Q, tauchen_trans_e, prev_ex_grid = res #unpack value function
    @unpack C, w, r, nϵ, nQ, R, Q_grid, ϵ_grid, n_periods, n_firms, ρ_q, σ_q, θ, α_n, n_sims = prim #unpack primitives

    firms_export_capital = ones(n_periods, n_firms, n_sims)
    firms_export_decisions = zeros(n_periods, n_firms, n_sims)
    firms_labor_decisions = ones(n_periods, n_firms, n_sims)
    firms_capital_decisions = ones(n_periods, n_firms, n_sims)
    firms_sales_non_exporter = zeros(n_periods, n_firms, n_sims)
    firms_sales = zeros(n_periods, n_firms, n_sims)
    firms_export_sales = zeros(n_periods, n_firms, n_sims)
    
    for k = 1:n_sims
        Random.seed!(k)
        for i = 2:n_periods

            Q[i,k] = exp(ρ_q*log(Q[i-1,k]) + rand(Normal(0, σ_q)))
            if Q[i,k] > 1.2
                Q[i,k] == 1.2
            elseif Q[i,k] < 0.8
                Q[i,k] == 0.8
            end
            Q_index = findmin(abs.(Q[i,k] .- Q_grid))[2]

            for j = 1:n_firms
                
                ϵ[i,j,k] = exp(ρ_e*log(ϵ[i-1,j,k]) + rand(Normal(0,σ_e)))
                ϵ_index = findmin(abs.(ϵ[i,j,k] .- ϵ_grid))[2]
                
                firms_export_decisions[i,j,k] = ex_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]

                if firms_export_decisions[i,j,k] == 0 && firms_export_capital[i-1,j,k] > 1
                    firms_export_capital[i,j,k] = firms_export_capital[i-1,j,k] - 1
                elseif firms_export_decisions[i,j,k] == 1
                    firms_export_capital[i,j,k] = n_prev_ex
                end
 
                firms_labor_decisions[i,j,k] = n_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]
                firms_capital_decisions[i,j,k] = k_func[Q_index, ϵ_index, floor(Int,firms_export_capital[i-1,j,k])]               
                firms_sales_non_exporter[i,j,k] = domestic_revenue(prim, res, [firms_labor_decisions[i,j,k] firms_capital_decisions[i,j,k] firms_export_decisions[i,j,k] ϵ_grid[ϵ_index] Q_grid[Q_index]])
                firms_sales[i,j,k] = total_revenue(prim, res, [firms_labor_decisions[i,j,k] firms_capital_decisions[i,j,k] firms_export_decisions[i,j,k] ϵ_grid[ϵ_index] Q_grid[Q_index]])
                firms_export_sales[i,j,k] = export_revenue(prim, res, [firms_labor_decisions[i,j,k] firms_capital_decisions[i,j,k] firms_export_decisions[i,j,k] ϵ_grid[ϵ_index] Q_grid[Q_index]])
            end
        end
    end

    return firms_export_decisions, firms_labor_decisions, firms_capital_decisions, firms_sales_non_exporter, firms_sales, firms_export_sales

end

function moment_calc(prim, res)
    firms_export_decisions, firms_labor_decisions, firms_capital_decisions, firms_sales_domestic, firms_sales = data_sim_delta(prim, res)
    val_annual_firms_export_decisions = zeros(12, prim.n_firms)
    annual_firms_export_decisions = zeros(12, prim.n_firms)
    annual_firms_export_change = zeros(12, prim.n_firms)
    annual_firms_sales_domestic = zeros(12, prim.n_firms)
    annual_firms_sales = zeros(12, prim.n_firms)

    for i = 1:48
        for j = 1:12
            if ceil(i/4) == j
                val_annual_firms_export_decisions[j,:] = val_annual_firms_export_decisions[j,:] .+ firms_export_decisions[i+400,:]
                annual_firms_sales_domestic[j,:] = annual_firms_sales_domestic[j,:] .+ firms_sales_domestic[i+400,:]
                annual_firms_sales[j,:] = annual_firms_sales[j,:] .+ firms_sales[i+400,:]
            end
        end
    end

    annual_firms_export_decisions = val_annual_firms_export_decisions .> 0

    for i = 2:12
        for j = 1:prim.n_firms
            annual_firms_export_change[i,j] = annual_firms_export_decisions[i,j]-annual_firms_export_decisions[i-1,j]
        end
    end

    annual_starter_rate = zeros(11)
    annual_stopper_rate = zeros(11)
    # Add 1 and 0 to export decisions once to make sure neither are ever 0
    annual_firms_export_decisions = hcat(annual_firms_export_decisions, ones(12), zeros(12))
    # Add 1 and -1 to export changes once to make sure neither are ever 0
    annual_firms_export_change = hcat(annual_firms_export_change, ones(12), -1 .* ones(12)) 
    exporter_i = countmap(annual_firms_export_change[1,:])
    exporter_i_minus1 = countmap(annual_firms_export_decisions[1,:])

    one_zero = AbstractVector{Int64}([1,0])
    one_negOne = AbstractVector{Int64}([1,-1])

    for i = 2:12
        exporter_i_minus1 = countmap(annual_firms_export_decisions[i-1,:])
        exporter_i = countmap(annual_firms_export_change[i,:])
        
        annual_starter_rate[i-1] = exporter_i[1]/exporter_i_minus1[0]
        annual_stopper_rate[i-1] = exporter_i[-1]/exporter_i_minus1[1]
    end

    sales_today = annual_firms_sales[2,:]
    sales_yesterday = annual_firms_sales[1,:]
    year = []
    firm_num = []
    reg_data = zeros(prim.n_firms, 3)
    for i = 3:12
        sales_today = [sales_today; annual_firms_sales[i,:]]
        sales_yesterday = [sales_yesterday; annual_firms_sales[i-1,:]]
    end

    for i = 2:12
        for j = 1:prim.n_firms
            year = [year; string("Y",i)]
            firm_num = [firm_num; string("F",j)]
        end
    end

    reg_DF = DataFrame(s_today = log.(sales_today), s_yesterday = log.(sales_yesterday), y = year, f = firm_num)
    reg_results = reg(reg_DF, @formula(s_today ~ s_yesterday + fe(y) + fe(f)))

    β_moment = coef(reg_results)

    output = zeros(5)

    output[1] = mean(annual_starter_rate)
    output[2] = mean(annual_stopper_rate)
    output[3] = mean(res.C_star ./ (res.Q[401:448,:] .^(-1*prim.θ)))
    output[4] = std(log.(annual_firms_sales_domestic[:,:]))/mean(log.(annual_firms_sales_domestic[:,:]))
    output[5] = β_moment[1]

    print(output)
    print("\n")
    print([prim.true_starter, prim.true_stopper, prim.true_ave_es_ratio, prim.true_coef_var, prim.true_a_exp_growth])
    print("\n")

    return error 
end

function moment_calc_nsims(firms_export_decisions, firms_labor_decisions, firms_capital_decisions, firms_sales_domestic, firms_sales)
    Random.seed!(17)
    val_annual_firms_export_decisions = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_export_decisions = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_export_change = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales_domestic = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_production = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales_foreign = zeros(12, prim.n_firms, prim.n_sims)
    annual_total_production = zeros(12)
    annual_total_sales_foreign = zeros(12)
    annual_total_sales = zeros(12)

    for j = 1:12
        val_annual_firms_export_decisions[j,:,:] = val_annual_firms_export_decisions[j,:,:] .+ firms_export_decisions[j+100,:,:]
        annual_firms_sales_domestic[j,:,:] = annual_firms_sales_domestic[j,:,:] .+ firms_sales_domestic[j+100,:,:]
        annual_firms_sales[j,:,:] = annual_firms_sales[j,:,:] .+ firms_sales[j+100,:,:]
        annual_firms_production[j,:,:] = annual_firms_production[j,:,:] .+ res.ϵ[j+100,:,:].*firms_labor_decisions[j+100,:,:].^prim.α_n .*firms_capital_decisions[j+100,:,:].^(1-prim.α_n)
        annual_firms_sales_foreign[j,:,:] = annual_firms_sales_foreign[j,:,:] .+ firms_sales[j+100,:,:] .- firms_sales_domestic[j+100,:,:]
    end
       
    annual_firms_export_decisions = val_annual_firms_export_decisions .> 0
    mean_annual_firms_sales_foreign = mean(annual_firms_sales_foreign, dims=3)
    mean_annual_firms_sales = mean(annual_firms_sales, dims=3)
    mean_annual_firms_production = mean(annual_firms_production, dims=3)

    for i = 1:12
        for j = 1:prim.n_firms
            if i == 1
                annual_total_sales_foreign[i] = annual_total_sales_foreign[i] + mean_annual_firms_sales_foreign[i,j]
                annual_total_sales[i] = annual_total_sales[i] + mean_annual_firms_sales[i,j]
                annual_total_production[i] = annual_total_production[i] + mean_annual_firms_production[i,j]
            else
                annual_total_sales_foreign[i] = annual_total_sales_foreign[i] + mean_annual_firms_sales_foreign[i,j] 
                annual_total_sales[i] = annual_total_sales[i] + mean_annual_firms_sales[i,j]
                annual_total_production[i] = annual_total_production[i] + mean_annual_firms_production[i,j]
                annual_firms_export_change[i,j,:] = annual_firms_export_decisions[i,j,:] .- annual_firms_export_decisions[i-1,j,:]
            end
        end
    end
    annual_starter_rate = zeros(11, prim.n_sims)
    annual_stopper_rate = zeros(11, prim.n_sims)
    # Add 1 and 0 to export decisions once to make sure neither are ever 0
    annual_firms_export_decisions = hcat(annual_firms_export_decisions, ones(12, 1, prim.n_sims), zeros(12, 1, prim.n_sims))
    # Add 1 and -1 to export changes once to make sure neither are ever 0
    annual_firms_export_change = hcat(annual_firms_export_change, ones(12, 1, prim.n_sims), -1 .* ones(12, 1, prim.n_sims)) 
    exporter_i = countmap(annual_firms_export_change[1,:,:])
    exporter_i_minus1 = countmap(annual_firms_export_decisions[1,:,:])

    for i = 2:12
        for k = 1:prim.n_sims
            exporter_i_minus1 = countmap(annual_firms_export_decisions[i-1,:,k])
            exporter_i = countmap(annual_firms_export_change[i,:,k])
            
            annual_starter_rate[i-1,k] = exporter_i[1]/exporter_i_minus1[0]
            annual_stopper_rate[i-1,k] = exporter_i[-1]/exporter_i_minus1[1]
        end
    end

    mean_annual_starter_rate = mean(annual_starter_rate, dims=2)
    mean_annual_stopper_rate = mean(annual_stopper_rate, dims=2)

    sales_today = mean(annual_firms_sales[2,:,:], dims=2)
    sales_yesterday = mean(annual_firms_sales[1,:,:], dims=2)
    year = []
    firm_num = []
    for i = 3:12
        sales_today = [sales_today; mean(annual_firms_sales[i,:,:],dims=2)]
        sales_yesterday = [sales_yesterday; mean(annual_firms_sales[i-1,:,:],dims=2)]
    end

    for i = 2:12
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

    output = zeros(5)

    output[1] = mean(annual_starter_rate)
    output[2] = mean(annual_stopper_rate)
    output[3] = mean(res.C_star ./ (res.Q[101:112,17] .^(-1*prim.θ)))
    output[4] = std(log.(annual_firms_sales_domestic[:,:,20]))/mean(log.(annual_firms_sales_domestic[:,:,20]))
    output[5] = β_moment[1]

    #print(output)
    #print("\n")
    #print([prim.true_starter, prim.true_stopper, prim.true_ave_es_ratio, prim.true_coef_var, prim.true_a_exp_growth])
    #print("\n")

    return mean_annual_starter_rate, mean_annual_stopper_rate, annual_total_sales_foreign, annual_total_sales, annual_total_production, output
end

function reentry_calcs(firms_export_decisions, firms_sales_domestic, firms_sales)
    val_annual_firms_export_decisions = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_export_decisions = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_export_change = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales_domestic = zeros(12, prim.n_firms, prim.n_sims)
    annual_firms_sales = zeros(12, prim.n_firms, prim.n_sims)
    prop_one_entry = zeros(prim.n_sims)
    prop_two_entry = zeros(prim.n_sims)
    prop_three_entry = zeros(prim.n_sims)
    prop_four_entry = zeros(prim.n_sims)
    prop_five_entry = zeros(prim.n_sims)
    prop_six_entry = zeros(prim.n_sims)
    prop_seven_entry = zeros(prim.n_sims)

    for j = 1:12
        val_annual_firms_export_decisions[j,:,:] = val_annual_firms_export_decisions[j,:,:] .+ firms_export_decisions[j+100,:,:]
        annual_firms_sales_domestic[j,:,:] = annual_firms_sales_domestic[j,:,:] .+ firms_sales_domestic[j+100,:,:]
        annual_firms_sales[j,:,:] = annual_firms_sales[j,:,:] .+ firms_sales[j+100,:,:]
    end

    annual_firms_export_decisions = val_annual_firms_export_decisions .> 0
    annual_firms_export_change[1,:,:] = copy(annual_firms_export_decisions[1,:,:])

    for i = 2:12
        for j = 1:prim.n_firms
            annual_firms_export_change[i,j,:] = annual_firms_export_decisions[i,j,:] .- annual_firms_export_decisions[i-1,j,:]
        end
    end

    for k = 1:prim.n_sims
        n_no_entry = 0
        n_one_entry = 0
        n_two_entries = 0
        n_three_entries = 0
        n_four_entries = 0
        n_five_entries = 0
        n_six_entries = 0
        n_seven_entries = 0
        export_entries = Vector{Int64}()
        export_entries_no_zero = Vector{Int64}()

        for i = 1:prim.n_firms
            exporter_i_map = countmap(vcat(annual_firms_export_change[:,i,k],1))
            if exporter_i_map[1] == 2
                n_one_entry = n_one_entry + 1
                export_entries = push!(export_entries, 1)
                export_entries_no_zero = push!(export_entries_no_zero, 1)
            elseif exporter_i_map[1] == 3
                n_two_entries = n_two_entries + 1
                export_entries =  push!(export_entries, 2)
                export_entries_no_zero = push!(export_entries_no_zero, 2)
            elseif exporter_i_map[1] == 4
                n_three_entries = n_three_entries + 1
                export_entries =  push!(export_entries, 3)
                export_entries_no_zero = push!(export_entries_no_zero, 3)
            elseif exporter_i_map[1] == 5
                n_four_entries = n_four_entries + 1
                export_entries =  push!(export_entries, 4)
                export_entries_no_zero = push!(export_entries_no_zero, 4)
            elseif exporter_i_map[1] == 6
                n_five_entries = n_five_entries + 1
                export_entries =  push!(export_entries, 5)
                export_entries_no_zero = push!(export_entries_no_zero, 5)
            elseif exporter_i_map[1] == 7
                n_six_entries = n_six_entries + 1
                export_entries =  push!(export_entries, 6)
                export_entries_no_zero = push!(export_entries_no_zero, 6)
            elseif exporter_i_map[1] == 8
                n_seven_entries = n_seven_entries + 1
                export_entries =  push!(export_entries, 7)
                export_entries_no_zero = push!(export_entries_no_zero, 7)
            elseif exporter_i_map[1] == 1
                n_no_entry = n_no_entry + 1
                export_entries =  push!(export_entries, 0)
            end
        end

        prop_two_entry[k] = n_two_entries/(prim.n_firms-n_no_entry)
        prop_three_entry[k] = n_three_entries/(prim.n_firms-n_no_entry)
        prop_four_entry[k] = n_four_entries/(prim.n_firms-n_no_entry)
        prop_five_entry[k] = n_five_entries/(prim.n_firms-n_no_entry)
        prop_six_entry[k] = n_six_entries/(prim.n_firms-n_no_entry)
        prop_seven_entry[k] = n_seven_entries/(prim.n_firms-n_no_entry)
        prop_one_entry[k] = n_one_entry/(prim.n_firms-n_no_entry)
    end

    prop_one_entry_mean = round(mean(prop_one_entry), digits = 4)
    prop_two_entry_mean = round(mean(prop_two_entry), digits = 4)
    prop_three_entry_mean = round(mean(prop_three_entry), digits = 4)
    prop_four_entry_mean = round(mean(prop_four_entry), digits = 4)
    prop_five_entry_mean = round(mean(prop_five_entry), digits = 4)

    return prop_one_entry_mean, prop_two_entry_mean, prop_three_entry_mean, prop_four_entry_mean, prop_five_entry_mean
end

#Bellman Operator
function Bellman_subsidy(prim::Primitives,res::Results,subsidy_rate::Float64)
    @unpack val_func, ex_func, n_func, k_func, val_func_subsidy, ex_func_subsidy, n_func_subsidy, k_func_subsidy, τ, C_star, σ_e, ρ_e, FC_0, FC_1, prev_ex_grid, n_prev_ex, tauchen_trans_Q, tauchen_trans_e = res #unpack value function
    @unpack Q_grid, Q_max, ϵ_grid, ϵ_max, R, r, w, σ_q, nQ, Q_grid, nϵ, ϵ_grid, ρ_q, α_n, θ, C, γ_0, γ_1 = prim #unpack model primitives
    v_next = zeros(nQ,nϵ,n_prev_ex) #next guess of value function to fill

    for Q_index = 1:nQ   
        Q = Q_grid[Q_index] #value of Q
        for ϵ_index = 1:nϵ
            ϵ = ϵ_grid[ϵ_index] #value of ϵ
            for prev_ex_index in 1:n_prev_ex
                prev_ex = prev_ex_grid[prev_ex_index]
                candidate_max = -1e10 #bad candidate max
                
                n_prev = 0
                k_prev = 0
                Π_prev = 0

                for ex_decision in 0:1 #loop over possible export decision
                    Π(x) = -1*((1+ex_decision*(Q*(1-res.τ))^θ*C_star/C)^(1/θ)*C^(1/θ)*ϵ*x[1]^(α_n*(θ-1)/θ)*x[2]^((1-α_n)*(θ-1)/θ)-w*x[1]-r*x[2])
                    x0 = [0.1, 0.1]
                    lower = [0.01, 0.01]
                    upper = [Inf, Inf]
                    inner_optimizer = GradientDescent()
                    optim_Π = optimize(Π, lower, upper, x0, Fminbox(inner_optimizer))

                    n_opt = Optim.minimizer(optim_Π)[1]
                    k_opt = Optim.minimizer(optim_Π)[2]
                    Π = -1*Optim.minimum(optim_Π)

                    # Should this be different ? Adjust based on export choice?
                    continuation_val = 0
                    for ϵp_index in 1:nϵ
                        for Qp_index in 1:nQ
                            # Continuation Val is assuming no more subsidy
                            continuation_val += R*tauchen_trans_e[ϵ_index,ϵp_index]*tauchen_trans_Q[Q_index,Qp_index]*val_func[Qp_index,ϵp_index,prev_ex_index]
                        end
                    end
                    total_Π = Π + continuation_val
                    if ex_decision == 0
                        Π_prev = total_Π
                        n_prev = n_opt
                        k_prev = k_opt
                    else
                        if prev_ex_index == n_prev_ex
                            total_Π += -1*FC_1
                        else
                            total_Π += -1*((1-prev_ex)*FC_0+FC_1)
                        end
                        
                        # Subsidy pays a fraction of FC_0 if they export, but only up to their sunk cost 
                        Π_subsidy = total_Π + min(subsidy_rate*FC_0, (1-prev_ex)*FC_0)

                        # If they export anyways:
                        if total_Π > Π_prev
                            ex_func_subsidy[Q_index, ϵ_index, prev_ex_index] = 1
                            n_func_subsidy[Q_index, ϵ_index, prev_ex_index] = n_opt
                            k_func_subsidy[Q_index, ϵ_index, prev_ex_index] = k_opt
                            # Continuation val is no subsidy but exporting
                            if total_Π > candidate_max
                                candidate_max = total_Π
                            end
                        # If they only export bc of the subsidy, continuation val is not exporting:
                        elseif Π_subsidy > Π_prev
                            ex_func_subsidy[Q_index, ϵ_index, prev_ex_index] = 1
                            n_func_subsidy[Q_index, ϵ_index, prev_ex_index] = n_opt
                            k_func_subsidy[Q_index, ϵ_index, prev_ex_index] = k_opt
                            if Π_prev > candidate_max
                                # Continuation val is no subsidy and not exporting
                                # Lol jk trying something out.
                                candidate_max = Π_prev
                            end
                        # If even the subsidy doesn't get them to export:
                        else
                            ex_func_subsidy[Q_index, ϵ_index, prev_ex_index] = 0
                            n_func_subsidy[Q_index, ϵ_index, prev_ex_index] = n_prev
                            k_func_subsidy[Q_index, ϵ_index, prev_ex_index] = k_prev
                            if Π_prev > candidate_max
                                # Continuation val is not exporting
                                candidate_max = Π_prev
                            end
                        end
                    end

                end
                v_next[Q_index, ϵ_index, prev_ex_index] = candidate_max #update value function
            end
        end
    end
    v_next #return next guess of value function
end

#Value function iteration
function V_iterate_subsidy(prim::Primitives, res::Results, subsidy_rate; tol::Float64 = 1e-3, err::Float64 = 1000.0)
    n = 0 #counter

    while err>tol #begin iteration
        v_next = Bellman_subsidy(prim, res, subsidy_rate) #spit out new vectors
        err = abs.(maximum(v_next.-res.val_func_subsidy))/abs(v_next[prim.nQ, prim.nϵ, 1]) #reset error level
        res.val_func_subsidy = v_next #update value function
        #println("Error:", err)
        n+=1
    end
    println("Value function converged in ", n, " iterations.")
end