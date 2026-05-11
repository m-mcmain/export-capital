clear all

import delimited Panel_Sim_delta.csv
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

egen first_year = min(year), by(firm)
gen in_first_year = first_year == year
egen years_exported = total(exported), by(firm)

gen changed_export = abs(exported-exported_prev)
gen changed_export_direction = exported-exported_prev

gen changed_export_enter = 1  if changed_export_direction > 0 & changed_export_direction != .
gen changed_export_exit = 1  if changed_export_direction < 0 & changed_export_direction != .
replace changed_export_enter = 1 if changed_export_enter == . & in_first_year*exported
replace changed_export_enter = 1 if exported_prev == . & exported
replace changed_export_enter = 0 if changed_export_enter == .
replace changed_export_exit = 0 if changed_export_exit == . 

egen total_stints = total(changed_export_enter), by(firm)
by firm: gen current_stint = sum(changed_export_enter)
replace current_stint = current_stint*exported

egen total_exits = total(changed_export_exit), by(firm)
by firm: gen current_exit = sum(changed_export_exit)
replace current_exit = current_exit*(1-exported)
gen reentry_occurs = (current_exit < total_stints)*(current_exit>0)


*****************************************************************
************ Length of Stints With Sim Data Delta ***************
*****************************************************************
egen max_year_current_stint = max(year), by(firm current_stint)
egen min_year_current_stint = min(year), by(firm current_stint)
by firm year current_stint, sort: gen length_current_stint = max_year_current_stint - min_year_current_stint + 1
replace length_current_stint = length_current_stint*exported

preserve

keep firm total_stints length_current_stint current_stint years_exported
drop if length_current_stint == 0
duplicates drop
tab length_current_stint if total_stints == 1
tab length_current_stint if total_stints == 1 & length_current_stint < 12
tab length_current_stint if total_stints > 1
tab length_current_stint if total_stints > 1 & current_stint == 1
tab length_current_stint if total_stints > 1 & current_stint > 1
// graph twoway histogram length_current_stint, frac by(total_stints)
// graph export .\images\export_stints_all_est.jpg, replace
// graph twoway histogram length_current_stint if total_stints < 4, frac by(total_stints)
// graph export .\images\export_stints_less4_est.jpg, replace

restore

***************************************************************
************ Survival Rates With Sim Data Delta ***************
***************************************************************
preserve

egen reentrant = max(reentry), by(firm)

gen fchanged_export_exit = f.changed_export_exit
gen f2changed_export_exit = f2.changed_export_exit
gen f3changed_export_exit = f3.changed_export_exit
gen f4changed_export_exit = f4.changed_export_exit
gen f5changed_export_exit = f5.changed_export_exit

keep if changed_export_enter & f5changed_export_exit != .

* Everyone
tab fchanged_export_exit if current_stint == 1
tab f2changed_export_exit if current_stint == 1 & !fchanged_export_exit
tab f3changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit
tab f4changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit
tab f5changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit

* Drop firms who export the whole time and not re-entrants
drop if years_exported == 12
tab fchanged_export_exit if current_stint == 1 & !reentry_occurs
tab f2changed_export_exit if current_stint == 1 & !reentry_occurs & !fchanged_export_exit
tab f3changed_export_exit if current_stint == 1 & !reentry_occurs & !fchanged_export_exit & !f2changed_export_exit
tab f4changed_export_exit if current_stint == 1 & !reentry_occurs & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit
tab f5changed_export_exit if current_stint == 1 & !reentry_occurs & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit

* Re-Entrants
tab fchanged_export_exit if current_stint > 1
tab f2changed_export_exit if current_stint > 1 & !fchanged_export_exit
tab f3changed_export_exit if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit
tab f4changed_export_exit if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit
tab f5changed_export_exit if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit

* Re-Entrants First Entry
tab fchanged_export_exit if current_stint == 1 & reentrant
tab f2changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit
tab f3changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit & !f2changed_export_exit
tab f4changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit
tab f5changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit

restore

********************************************************
* Export Intensity, First-Time Entrants vs. Re-Entrants:
********************************************************
gen esr = export_sales/sales
gen log_esr = log(export_sales/sales)
gen reentry = current_stint > 1
gen permanent_exporter = years_exported == 12
egen reentrant = max(reentry), by(firm)
bysort firm current_stint: gen years_exporting = sum(exported)

tab years_exporting reentry if years_exported != 12, sum(esr) nost
tab years_exporting reentry if reentrant, sum(esr) nost
tab years_exporting if !reentrant & years_exported != 12, sum(esr) nost

reghdfe log_esr reentry permanent_exporter if changed_export_est_enter, absorb(ciiu4_max)
reghdfe log_esr reentry permanent_exporter paper_ind food_ind textile_ind if changed_export_est_enter & R_T
reghdfe log_esr reentry if changed_export_est_enter & reentrant, absorb(ciiu4_max)