clear all

********************************************************
************ Redprob With Sim Data Delta ***************
********************************************************

import delimited "../../export-capital sim data/Panel_Sim_delta.csv"
rename (v1 v2 v3 v4 v5 v6) (firm year exported capital sales export_sales)
gen lsales = log(sales)
gen lcapital = log(capital)

xtset firm year
gen lsales_lag = l.lsales
gen lcapital_lag = l.lcapital

gen exported_prev = l.exported
gen exported_prev2 = l2.exported*(1-exported_prev)
gen exported_prev3 = l3.exported*(1-exported_prev2)*(1-exported_prev)
gen exported_prev4 = l4.exported*(1-exported_prev3)*(1-exported_prev2)*(1-exported_prev)

*** Year Dummies:
forval i=1(1)12{
	gen y_`i' = year==`i'
}
gen year_n_5 = year-4

*** Var Macros
global lags "exported_prev exported_prev2 exported_prev3 exported_prev4"
* exported_est_prev6"
preserve

keep if _n < 2827

xtset firm year_n_5

redprob exported $lags lsales_lag y_6 y_7 y_8 y_9 y_10 y_11 y_12 (lsales_lag) if year >= 5, i(firm) t(year_n_5) quad(24)

restore

***********************************************************
************ Redprob With Sim Data no Delta ***************
***********************************************************
clear all
import delimited "../../export-capital sim data/Panel_Sim_noDelta.csv"
rename (v1 v2 v3 v4 v5) (firm year exported capital sales)
gen lsales = log(sales)
gen lcapital = log(capital)

xtset firm year
gen lsales_lag = l.lsales
gen lcapital_lag = l.lcapital

gen exported_prev = l.exported
gen exported_prev2 = l2.exported*(1-exported_prev)
gen exported_prev3 = l3.exported*(1-exported_prev2)*(1-exported_prev)
gen exported_prev4 = l4.exported*(1-exported_prev3)*(1-exported_prev2)*(1-exported_prev)

*** Year Dummies:
forval i=1(1)12{
	gen y_`i' = year==`i'
}
gen year_n_5 = year-4

*** Var Macros
global lags "exported_prev exported_prev2 exported_prev3 exported_prev4"
* exported_est_prev6"
preserve

*keep if _n < 339121

xtset firm year_n_5

redprob exported $lags lsales_lag y_6 y_7 y_8 y_9 y_10 y_11 y_12 (lsales_lag) if year >= 5, i(firm) t(year_n_5) quadrat(24)

restore