clear all

import delimited "../../export-capital sim data/Panel_Sim_delta.csv"
rename (v1 v2 v3 v4 v5 v6 v7) (firm year exported capital sales export_sales productivity)
gen lsales = log(sales)
gen lcapital = log(capital)

* Drop Duplicates for Now
duplicates drop firm year, force

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
tab length_current_stint if total_stints == 1, matcell(frequencies)
matrix define length_stints_all_delta = frequencies/r(N)

tab length_current_stint if total_stints == 1 & length_current_stint < 12, matcell(frequencies)
matrix define length_stints_noForever_delta = frequencies/r(N)

tab length_current_stint if total_stints > 1

tab length_current_stint if total_stints > 1 & current_stint == 1, matcell(frequencies)
matrix define length_stints_reentry_first_delt = frequencies/r(N)

tab length_current_stint if total_stints > 1 & current_stint > 1, matcell(frequencies)
matrix define length_stints_reentries_delta = frequencies/r(N)

matrix coljoinbyname length_stints_delta = length_stints_all_delta length_stints_noForever_delta length_stints_reentry_first_delt length_stints_reentries_delta
// // graph twoway histogram length_current_stint, frac by(total_stints)
// // graph export .\images\export_stints_all_est.jpg, replace
// // graph twoway histogram length_current_stint if total_stints < 4, frac by(total_stints)
// // graph export .\images\export_stints_less4_est.jpg, replace

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

* No Re-Entrants
matrix define survival_all = (0 \ 0 \ 0 \ 0 \ 0)
matrix colnames survival_all = "One-Time Entry"

tab fchanged_export_exit if current_stint == 1 & !reentry_occurs, matcell(frequencies)
matrix survival_all[1,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f2changed_export_exit if current_stint == 1 & !fchanged_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_all[2,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f3changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_all[3,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f4changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_all[4,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f5changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_all[5,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

* Drop firms who export the whole time and not re-entrants
drop if years_exported == 12
matrix define survival_noForever = (0 \ 0 \ 0 \ 0 \ 0)
matrix colnames survival_noForever = "One-Time Entry No Permanent"

tab fchanged_export_exit if current_stint == 1 & !reentry_occurs, matcell(frequencies)
matrix survival_noForever[1,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f2changed_export_exit if current_stint == 1 & !fchanged_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_noForever[2,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f3changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_noForever[3,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f4changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_noForever[4,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f5changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_noForever[5,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

* Re-Entrants
matrix define survival_reentry = (0 \ 0 \ 0 \ 0 \ 0)
matrix colnames survival_reentry = "Re-Entry"

tab fchanged_export_exit if current_stint > 1, matcell(frequencies)
matrix survival_reentry [1,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f2changed_export_exit if current_stint > 1 & !fchanged_export_exit, matcell(frequencies)
matrix survival_reentry[2,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f3changed_export_exit if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit, matcell(frequencies)
matrix survival_reentry[3,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f4changed_export_exit if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit, matcell(frequencies)
matrix survival_reentry[4,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f5changed_export_exit if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit, matcell(frequencies)
matrix survival_reentry[5,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

* Re-Entrants Productivity
sum productivity if !fchanged_export_exit & current_stint > 1, det

sum productivity if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit, det

sum productivity if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit, det

sum productivity if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit, det

sum productivity if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit & !f5changed_export_exit, det

* Re-Entrants First Entry
matrix define survival_reentry_first = (0 \ 0 \ 0 \ 0 \ 0)
matrix colnames survival_reentry_first = "Re-Entry First"

tab fchanged_export_exit if current_stint == 1 & reentrant
matrix survival_reentry_first[1,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f2changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit, matcell(frequencies)
matrix survival_reentry_first[2,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f3changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit & !f2changed_export_exit, matcell(frequencies)
matrix survival_reentry_first[3,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f4changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit, matcell(frequencies)
matrix survival_reentry_first[4,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f5changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit, matcell(frequencies)
matrix survival_reentry_first[5,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

matrix coljoinbyname survival_delta = survival_all survival_noForever survival_reentry_first survival_reentry

restore

**************************************************************
* Export Intensity, First-Time Entrants vs. Re-Entrants Delta:
**************************************************************
gen esr = export_sales/sales
gen log_esr = log(export_sales/sales)
gen reentry = current_stint > 1
gen permanent_exporter = years_exported == 12
egen reentrant = max(reentry), by(firm)
bysort firm current_stint: gen years_exporting = sum(exported)

// tab years_exporting reentry if years_exported != 12, sum(esr) nost
// tab years_exporting reentry if reentrant, sum(esr) nost
// tab years_exporting if !reentrant & years_exported != 12, sum(esr) nost

* Total:
tab years_exporting, sum(esr)

* Split by Exporter Type
tabstat esr if !reentry_occurs & length_current_stint == 1, by(years_exporting) save
matrix rowjoinbyname esr_y1 = r(Stat1)
tabstat esr if !reentry_occurs & length_current_stint == 2, by(years_exporting) save
matrix rowjoinbyname esr_y2 = r(Stat1) r(Stat2)
tabstat esr if !reentry_occurs & length_current_stint == 3, by(years_exporting) save
matrix rowjoinbyname esr_y3 = r(Stat1) r(Stat2) r(Stat3)
tabstat esr if !reentry_occurs & length_current_stint == 4, by(years_exporting) save
matrix rowjoinbyname esr_y4 = r(Stat1) r(Stat2) r(Stat3) r(Stat4)
tabstat esr if !reentry_occurs & length_current_stint == 5, by(years_exporting) save
matrix rowjoinbyname esr_y5 = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
tabstat esr if !reentry_occurs & length_current_stint > 5, by(years_exporting) save
matrix rowjoinbyname esr_y5plus = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
matrix coljoinbyname esr_one = esr_y1 esr_y2 esr_y3 esr_y4 esr_y5 esr_y5plus
* Plot!
preserve

clear
svmat esr_one

label variable esr_one1 "1 yr spell"
label variable esr_one2 "2"
label variable esr_one3 "3"
label variable esr_one4 "4"
label variable esr_one5 "5"
label variable esr_one6 "5+"
gen years_exporting = _n
label variable years_exporting "Years Exporting"

graph twoway (line (esr_one1 esr_one2 esr_one3 esr_one4 esr_one5 esr_one6) years_exporting, ytitle("Export to Sales Ratio") legend(position(0) bplacement(seast) order(1 2 3 4 5 6)) lwidth(thick thick thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "148 203 236" "220 205 125" "221 221 221") msymbol(O)) (scatter esr_one1 years_exporting, mcolor("51 34 136"))
graph export ./images/esr/esr_one_delta.pdf, replace

restore

tabstat esr if !reentry_occurs & length_current_stint == 1 & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y1 = r(Stat1) 
tabstat esr if !reentry_occurs & length_current_stint == 2 & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y2 = r(Stat1) r(Stat2)
tabstat esr if !reentry_occurs & length_current_stint == 3 & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y3 = r(Stat1) r(Stat2) r(Stat3)
tabstat esr if !reentry_occurs & length_current_stint == 4 & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y4 = r(Stat1) r(Stat2) r(Stat3) r(Stat4)
tabstat esr if !reentry_occurs & length_current_stint == 5 & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y5 = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
tabstat esr if !reentry_occurs & length_current_stint >5  & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y5plus = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
matrix coljoinbyname esr_one_noperm = esr_y1 esr_y2 esr_y3 esr_y4 esr_y5 esr_y5plus
* Plot!
preserve

clear
svmat esr_one_noperm

label variable esr_one_noperm1 "1 yr spell"
label variable esr_one_noperm2 "2"
label variable esr_one_noperm3 "3"
label variable esr_one_noperm4 "4"
label variable esr_one_noperm5 "5"
label variable esr_one_noperm6 "5+"
gen years_exporting = _n
label variable years_exporting "Years Exporting"

graph twoway (line (esr_one_noperm1 esr_one_noperm2 esr_one_noperm3 esr_one_noperm4 esr_one_noperm5 esr_one_noperm6) years_exporting, ytitle("Export to Sales Ratio") legend(position(0) bplacement(seast) order(1 2 3 4 5 6)) lwidth(thick thick thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "148 203 236" "220 205 125" "221 221 221") msymbol(O)) (scatter esr_one_noperm1 years_exporting, mcolor("51 34 136"))
graph export ./images/esr/esr_one_noperm_delta.pdf, replace

restore

tabstat esr if current_stint == 1 & reentrant & length_current_stint == 1, by(years_exporting) save
matrix rowjoinbyname esr_y1 = r(Stat1)
tabstat esr if current_stint == 1 & reentrant & length_current_stint == 2, by(years_exporting) save
matrix rowjoinbyname esr_y2 = r(Stat1) r(Stat2)
tabstat esr if current_stint == 1 & reentrant & length_current_stint == 3, by(years_exporting) save
matrix rowjoinbyname esr_y3 = r(Stat1) r(Stat2) r(Stat3)
tabstat esr if current_stint == 1 & reentrant & length_current_stint == 4, by(years_exporting) save
matrix rowjoinbyname esr_y4 = r(Stat1) r(Stat2) r(Stat3) r(Stat4)
tabstat esr if current_stint == 1 & reentrant & length_current_stint == 5, by(years_exporting) save
matrix rowjoinbyname esr_y5 = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
tabstat esr if current_stint == 1 & reentrant & length_current_stint > 5, by(years_exporting) save
matrix rowjoinbyname esr_y5plus = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
matrix coljoinbyname esr_reentrant_first = esr_y1 esr_y2 esr_y3 esr_y4 esr_y5 esr_y5plus
* Plot!
preserve

clear
svmat esr_reentrant_first

label variable esr_reentrant_first1 "1 yr spell"
label variable esr_reentrant_first2 "2"
label variable esr_reentrant_first3 "3"
label variable esr_reentrant_first4 "4"
label variable esr_reentrant_first5 "5"
label variable esr_reentrant_first6 "5+"
gen years_exporting = _n
label variable years_exporting "Years Exporting"

graph twoway (line (esr_reentrant_first1 esr_reentrant_first2 esr_reentrant_first3 esr_reentrant_first4 esr_reentrant_first5 esr_reentrant_first6) years_exporting, ytitle("Export to Sales Ratio") legend(position(0) bplacement(seast) order(1 2 3 4 5 6)) lwidth(thick thick thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "148 203 236" "220 205 125" "221 221 221") msymbol(O)) (scatter esr_reentrant_first1 years_exporting, mcolor("51 34 136"))
graph export ./images/esr/esr_reentrant_first_delta.pdf, replace

restore


tabstat esr if current_stint > 1 & reentrant & length_current_stint == 1, by(years_exporting) save
matrix rowjoinbyname esr_y1 = r(Stat1)
tabstat esr if current_stint > 1 & reentrant & length_current_stint == 2, by(years_exporting) save
matrix rowjoinbyname esr_y2 = r(Stat1) r(Stat2)
tabstat esr if current_stint > 1 & reentrant & length_current_stint == 3, by(years_exporting) save
matrix rowjoinbyname esr_y3 = r(Stat1) r(Stat2) r(Stat3)
tabstat esr if current_stint > 1 & reentrant & length_current_stint == 4, by(years_exporting) save
matrix rowjoinbyname esr_y4 = r(Stat1) r(Stat2) r(Stat3) r(Stat4)
tabstat esr if current_stint > 1 & reentrant & length_current_stint == 5, by(years_exporting) save
matrix rowjoinbyname esr_y5 = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
tabstat esr if current_stint > 1 & reentrant & length_current_stint > 5, by(years_exporting) save
matrix rowjoinbyname esr_y5plus = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
matrix coljoinbyname esr_reentries = esr_y1 esr_y2 esr_y3 esr_y4 esr_y5 esr_y5plus
* Plot!
preserve

clear
svmat esr_reentries

label variable esr_reentries1 "1 yr spell"
label variable esr_reentries2 "2"
label variable esr_reentries3 "3"
label variable esr_reentries4 "4"
label variable esr_reentries5 "5"
label variable esr_reentries6 "5+"
gen years_exporting = _n
label variable years_exporting "Years Exporting"

graph twoway (line (esr_reentries1 esr_reentries2 esr_reentries3 esr_reentries4 esr_reentries5 esr_reentries6) years_exporting, ytitle("Export to Sales Ratio") legend(position(0) bplacement(seast) order(1 2 3 4 5 6)) lwidth(thick thick thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "148 203 236" "220 205 125" "221 221 221")) (scatter esr_reentries1 years_exporting, mcolor("51 34 136"))
graph export ./images/esr/esr_reentries_delta.pdf, replace

restore

// reghdfe log_esr reentry permanent_exporter if changed_export_est_enter, absorb(ciiu4_max)
// reghdfe log_esr reentry permanent_exporter paper_ind food_ind textile_ind if changed_export_est_enter & R_T
// reghdfe log_esr reentry if changed_export_est_enter & reentrant, absorb(ciiu4_max)


********************************************************************************
********************************************************************************
* 					Same thing but for Sunk Cost Data						   *
********************************************************************************
********************************************************************************
clear

import delimited "../../export-capital sim data/Panel_Sim_noDelta.csv"
rename (v1 v2 v3 v4 v5 v6 v7) (firm year exported capital sales export_sales productivity)
gen lsales = log(sales)
gen lcapital = log(capital)

* Drop Duplicates for Now
duplicates drop firm year, force

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
********** Length of Stints With Sim Data No Delta **************
*****************************************************************
egen max_year_current_stint = max(year), by(firm current_stint)
egen min_year_current_stint = min(year), by(firm current_stint)
by firm year current_stint, sort: gen length_current_stint = max_year_current_stint - min_year_current_stint + 1
replace length_current_stint = length_current_stint*exported

preserve

keep firm total_stints length_current_stint current_stint years_exported
drop if length_current_stint == 0
duplicates drop
tab length_current_stint if total_stints == 1, matcell(frequencies)
matrix define length_stints_all_nodelta = frequencies/r(N)

tab length_current_stint if total_stints == 1 & length_current_stint < 12, matcell(frequencies)
matrix define length_stints_noForever_nodelta = frequencies/r(N)

tab length_current_stint if total_stints > 1

tab length_current_stint if total_stints > 1 & current_stint == 1, matcell(frequencies)
matrix define length_stints_reentry_first_node = frequencies/r(N)

tab length_current_stint if total_stints > 1 & current_stint > 1, matcell(frequencies)
matrix define length_stints_reentries_nodelta = frequencies/r(N)

matrix coljoinbyname length_stints_nodelta = length_stints_all_nodelta length_stints_noForever_nodelta length_stints_reentry_first_node length_stints_reentries_nodelta
// graph twoway histogram length_current_stint, frac by(total_stints)
// graph export .\images\export_stints_all_est.jpg, replace
// graph twoway histogram length_current_stint if total_stints < 4, frac by(total_stints)
// graph export .\images\export_stints_less4_est.jpg, replace

restore

***************************************************************
********** Survival Rates With Sim Data No Delta **************
***************************************************************
preserve

egen reentrant = max(reentry), by(firm)

gen fchanged_export_exit = f.changed_export_exit
gen f2changed_export_exit = f2.changed_export_exit
gen f3changed_export_exit = f3.changed_export_exit
gen f4changed_export_exit = f4.changed_export_exit
gen f5changed_export_exit = f5.changed_export_exit

keep if changed_export_enter & f5changed_export_exit != .

* No Re-Entrants
matrix define survival_all = (0 \ 0 \ 0 \ 0 \ 0)
matrix colnames survival_all = "One-Time Entry"

tab fchanged_export_exit if current_stint == 1 & !reentry_occurs, matcell(frequencies)
matrix survival_all[1,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f2changed_export_exit if current_stint == 1 & !fchanged_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_all[2,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f3changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_all[3,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f4changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_all[4,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f5changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_all[5,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

* Drop firms who export the whole time and not re-entrants
drop if years_exported == 12
matrix define survival_noForever = (0 \ 0 \ 0 \ 0 \ 0)
matrix colnames survival_noForever = "One-Time Entry No Permanent"

tab fchanged_export_exit if current_stint == 1 & !reentry_occurs, matcell(frequencies)
matrix survival_noForever[1,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f2changed_export_exit if current_stint == 1 & !fchanged_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_noForever[2,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f3changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_noForever[3,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f4changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_noForever[4,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f5changed_export_exit if current_stint == 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit & !reentry_occurs, matcell(frequencies)
matrix survival_noForever[5,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

* Re-Entrants
matrix define survival_reentry = (0 \ 0 \ 0 \ 0 \ 0)
matrix colnames survival_reentry = "Re-Entry"

tab fchanged_export_exit if current_stint > 1, matcell(frequencies)
matrix survival_reentry[1,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])
if frequencies[2,1] == .{
	matrix survival_reentry[1,1] = 1 
}

tab f2changed_export_exit if current_stint > 1 & !fchanged_export_exit, matcell(frequencies)
matrix survival_reentry[2,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f3changed_export_exit if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit, matcell(frequencies)
matrix survival_reentry[3,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f4changed_export_exit if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit, matcell(frequencies)
matrix survival_reentry[4,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f5changed_export_exit if current_stint > 1 & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit, matcell(frequencies)
matrix survival_reentry[5,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

* Re-Entrants First Entry
matrix define survival_reentry_first = (0 \ 0 \ 0 \ 0 \ 0)
matrix colnames survival_reentry_first = "Re-Entry First"

tab fchanged_export_exit if current_stint == 1 & reentrant
matrix survival_reentry_first[1,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f2changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit, matcell(frequencies)
matrix survival_reentry_first[2,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f3changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit & !f2changed_export_exit, matcell(frequencies)
matrix survival_reentry_first[3,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f4changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit, matcell(frequencies)
matrix survival_reentry_first[4,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

tab f5changed_export_exit if current_stint == 1 & reentrant & !fchanged_export_exit & !f2changed_export_exit & !f3changed_export_exit & !f4changed_export_exit, matcell(frequencies)
matrix survival_reentry_first[5,1] = frequencies[1,1]/(frequencies[1,1]+frequencies[2,1])

matrix coljoinbyname survival_nodelta = survival_all survival_noForever survival_reentry_first survival_reentry

restore

**************************************************************
* Export Intensity, First-Time Entrants vs. Re-Entrants Delta:
**************************************************************
gen esr = export_sales/sales
gen log_esr = log(export_sales/sales)
gen reentry = current_stint > 1
gen permanent_exporter = years_exported == 12
egen reentrant = max(reentry), by(firm)
bysort firm current_stint: gen years_exporting = sum(exported)

// tab years_exporting reentry if years_exported != 12, sum(esr) nost
// tab years_exporting reentry if reentrant, sum(esr) nost
// tab years_exporting if !reentrant & years_exported != 12, sum(esr) nost

* Total:
tab years_exporting, sum(esr)

tabstat esr if !reentry_occurs & length_current_stint == 1, by(years_exporting) save
matrix rowjoinbyname esr_y1 = r(Stat1)
tabstat esr if !reentry_occurs & length_current_stint == 2, by(years_exporting) save
matrix rowjoinbyname esr_y2 = r(Stat1) r(Stat2)
tabstat esr if !reentry_occurs & length_current_stint == 3, by(years_exporting) save
matrix rowjoinbyname esr_y3 = r(Stat1) r(Stat2) r(Stat3)
tabstat esr if !reentry_occurs & length_current_stint == 4, by(years_exporting) save
matrix rowjoinbyname esr_y4 = r(Stat1) r(Stat2) r(Stat3) r(Stat4)
tabstat esr if !reentry_occurs & length_current_stint == 5, by(years_exporting) save
matrix rowjoinbyname esr_y5 = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
tabstat esr if !reentry_occurs & length_current_stint > 5, by(years_exporting) save
matrix rowjoinbyname esr_y5plus = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
matrix coljoinbyname esr_one = esr_y1 esr_y2 esr_y3 esr_y4 esr_y5 esr_y5plus
* Plot!
preserve

clear
svmat esr_one

label variable esr_one1 "1 yr spell"
label variable esr_one2 "2"
label variable esr_one3 "3"
label variable esr_one4 "4"
label variable esr_one5 "5"
label variable esr_one6 "5+"
gen years_exporting = _n
label variable years_exporting "Years Exporting"

graph twoway (line (esr_one1 esr_one2 esr_one3 esr_one4 esr_one5 esr_one6) years_exporting, ytitle("Export to Sales Ratio") legend(position(0) bplacement(seast) order(1 2 3 4 5 6)) lwidth(thick thick thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "148 203 236" "220 205 125" "221 221 221") msymbol(O)) (scatter esr_one1 years_exporting, mcolor("51 34 136"))
graph export ./images/esr/esr_one_nodelta.pdf, replace

restore

tabstat esr if !reentry_occurs & length_current_stint == 1 & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y1 = r(Stat1) 
tabstat esr if !reentry_occurs & length_current_stint == 2 & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y2 = r(Stat1) r(Stat2)
tabstat esr if !reentry_occurs & length_current_stint == 3 & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y3 = r(Stat1) r(Stat2) r(Stat3)
tabstat esr if !reentry_occurs & length_current_stint == 4 & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y4 = r(Stat1) r(Stat2) r(Stat3) r(Stat4)
tabstat esr if !reentry_occurs & length_current_stint == 5 & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y5 = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
tabstat esr if !reentry_occurs & length_current_stint >5  & years_exported != 12, by(years_exporting) save
matrix rowjoinbyname esr_y5plus = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
matrix coljoinbyname esr_one_noperm = esr_y1 esr_y2 esr_y3 esr_y4 esr_y5 esr_y5plus
* Plot!
preserve

clear
svmat esr_one_noperm

label variable esr_one_noperm1 "1 yr spell"
label variable esr_one_noperm2 "2"
label variable esr_one_noperm3 "3"
label variable esr_one_noperm4 "4"
label variable esr_one_noperm5 "5"
label variable esr_one_noperm6 "5+"
gen years_exporting = _n
label variable years_exporting "Years Exporting"

graph twoway (line (esr_one_noperm1 esr_one_noperm2 esr_one_noperm3 esr_one_noperm4 esr_one_noperm5 esr_one_noperm6) years_exporting, ytitle("Export to Sales Ratio") legend(position(0) bplacement(seast) order(1 2 3 4 5 6)) lwidth(thick thick thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "148 203 236" "220 205 125" "221 221 221") msymbol(O)) (scatter esr_one_noperm1 years_exporting, mcolor("51 34 136"))
graph export ./images/esr/esr_one_noperm_nodelta.pdf, replace

restore


tabstat esr if current_stint == 1 & reentrant & length_current_stint == 1, by(years_exporting) save
matrix rowjoinbyname esr_y1 = r(Stat1)
tabstat esr if current_stint == 1 & reentrant & length_current_stint == 2, by(years_exporting) save
matrix rowjoinbyname esr_y2 = r(Stat1) r(Stat2)
tabstat esr if current_stint == 1 & reentrant & length_current_stint == 3, by(years_exporting) save
matrix rowjoinbyname esr_y3 = r(Stat1) r(Stat2) r(Stat3)
tabstat esr if current_stint == 1 & reentrant & length_current_stint == 4, by(years_exporting) save
matrix rowjoinbyname esr_y4 = r(Stat1) r(Stat2) r(Stat3) r(Stat4)
tabstat esr if current_stint == 1 & reentrant & length_current_stint == 5, by(years_exporting) save
matrix rowjoinbyname esr_y5 = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
tabstat esr if current_stint == 1 & reentrant & length_current_stint > 5, by(years_exporting) save
matrix rowjoinbyname esr_y5plus = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
matrix coljoinbyname esr_reentrant_first = esr_y1 esr_y2 esr_y3 esr_y4 esr_y5 esr_y5plus
* Plot!
preserve

clear
svmat esr_reentrant_first

label variable esr_reentrant_first1 "1 yr spell"
label variable esr_reentrant_first2 "2"
label variable esr_reentrant_first3 "3"
label variable esr_reentrant_first4 "4"
label variable esr_reentrant_first5 "5"
label variable esr_reentrant_first6 "5+"
gen years_exporting = _n
label variable years_exporting "Years Exporting"

graph twoway (line (esr_reentrant_first1 esr_reentrant_first2 esr_reentrant_first3 esr_reentrant_first4 esr_reentrant_first5 esr_reentrant_first6) years_exporting, ytitle("Export to Sales Ratio") legend(position(0) bplacement(seast) order(1 2 3 4 5 6)) lwidth(thick thick thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "148 203 236" "220 205 125" "221 221 221") msymbol(O)) (scatter esr_reentrant_first1 years_exporting, mcolor("51 34 136"))
graph export ./images/esr/esr_reentrant_first_nodelta.pdf, replace

restore


tabstat esr if current_stint > 1 & reentrant & length_current_stint == 1, by(years_exporting) save
matrix rowjoinbyname esr_y1 = r(Stat1)
tabstat esr if current_stint > 1 & reentrant & length_current_stint == 2, by(years_exporting) save
matrix rowjoinbyname esr_y2 = r(Stat1) r(Stat2)
tabstat esr if current_stint > 1 & reentrant & length_current_stint == 3, by(years_exporting) save
matrix rowjoinbyname esr_y3 = r(Stat1) r(Stat2) r(Stat3)
tabstat esr if current_stint > 1 & reentrant & length_current_stint == 4, by(years_exporting) save
matrix rowjoinbyname esr_y4 = r(Stat1) r(Stat2) r(Stat3) r(Stat4)
tabstat esr if current_stint > 1 & reentrant & length_current_stint == 5, by(years_exporting) save
matrix rowjoinbyname esr_y5 = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
tabstat esr if current_stint > 1 & reentrant & length_current_stint > 5, by(years_exporting) save
matrix rowjoinbyname esr_y5plus = r(Stat1) r(Stat2) r(Stat3) r(Stat4) r(Stat5)
matrix coljoinbyname esr_reentries = esr_y1 esr_y2 esr_y3 esr_y4 esr_y5 esr_y5plus
* Plot!
preserve

clear
svmat esr_reentries

label variable esr_reentries1 "1 yr spell"
label variable esr_reentries2 "2"
label variable esr_reentries3 "3"
label variable esr_reentries4 "4"
label variable esr_reentries5 "5"
label variable esr_reentries6 "5+"
gen years_exporting = _n
label variable years_exporting "Years Exporting"

graph twoway (line (esr_reentries1 esr_reentries2 esr_reentries3 esr_reentries4 esr_reentries5 esr_reentries6) years_exporting, ytitle("Export to Sales Ratio") legend(position(0) bplacement(seast) order(1 2 3 4 5 6)) lwidth(thick thick thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "148 203 236" "220 205 125" "221 221 221")) (scatter esr_reentries1 years_exporting, mcolor("51 34 136"))
graph export ./images/esr/esr_reentries_nodelta.pdf, replace

restore


********************************************************************************
********************************************************************************
* 					            Creat Plots       						       *
********************************************************************************
********************************************************************************

**************************************
*           Length Stints            *
**************************************
clear
svmat length_stints_delta

label variable length_stints_delta1 "Export Capital"
label variable length_stints_delta2 "Export Capital"
label variable length_stints_delta3 "Export Capital"
label variable length_stints_delta4 "Export Capital"

svmat length_stints_nodelta

label variable length_stints_nodelta1 "Sunk Cost"
label variable length_stints_nodelta2 "Sunk Cost"
label variable length_stints_nodelta3 "Sunk Cost"
label variable length_stints_nodelta4 "Sunk Cost"

matrix define length_stints_data = (0.159, 0.383, 0.426, 0.375 \ 0.058, 0.139, 0.177, 0.228 \ 0.03, 0.073, 0.112, 0.084 \ 0.033, 0.079, 0.089, 0.086 \ 0.017, 0.041, 0.094, 0.107 \ 0.005, 0.013, 0.0065, 0.046 \ 0.011, 0.026, 0.046, 0.027 \ 0.011, 0.026, 0.029, 0.005 \ 0.017, 0.041, 0.018, 0.025 \ 0.029, 0.069, 0.013, 0.018 \ 0.046, 0.111, ., . \ 0.585, ., ., .)
svmat length_stints_data

label variable length_stints_data1 "One-Time Entrants"
label variable length_stints_data2 "No Permanent One-Time Entrants"
label variable length_stints_data3 "Re-Entrants, First"
label variable length_stints_data4 "Re-Entries"

matrix define length_stints_post_y1_data = (0.395, 0.588, 0.375 \ 0.134, 0.226, 0.228 \ 0.079, 0.092, 0.084 \ 0.105, 0.039, 0.086 \ 0.054, 0.016, 0.107 \ 0.017, 0.016, 0.046 \ 0.003, 0.013, 0.027 \ 0.009, 0.007, 0.005 \ 0.04, 0.003, 0.025 \ 0.063, ., 0.018 \ 0.108, ., .)
svmat length_stints_post_y1_data

label variable length_stints_post_y1_data1 "One-Time Entrants"
label variable length_stints_post_y1_data2 "Re-Entrants, First"
label variable length_stints_post_y1_data3 "Re-Entries"

gen years = _n
**********************
* Length Stints Data *
**********************
reshape long length_stints_data, i(years) j(data_type)
label define length_stints_data_lab 1 "One-Time Entrants" 2 "No Permanent One-Time Entrants" 3 "Re-Entrants First" 4 "Re-Entries"
label values data_type length_stints_data_lab
graph bar length_stints_data, over(data_type) over(years) asyvars ytitle("Percent") b1title("Years") legend(position(0) bplacement(nwest)) bar(1, color(51 34 136)) bar(2, color(17 119 51)) bar(3, color(204 102 119)) bar(4, color(82 137 199)) ysc(r(0 0.6))
graph export ./images/length_stints_data.pdf, replace as(pdf)

reshape wide
******************************
* Length Stints Post Y1 Data *
******************************
reshape long length_stints_post_y1_data, i(years) j(data_type)
*label define length_stints_post_y1_data_lab 1 "One-Time Entrants" 2 "Re-Entrants First" 3 "Re-Entries"
label values data_type length_stints_post_y1_data_lab
graph bar length_stints_post_y1_data, over(data_type) over(years) asyvars ytitle("Percent") b1title("Years") legend(position(0) bplacement(neast)) bar(1, color(17 119 51)) bar(2, color(204 102 119)) bar(3, color(82 137 199)) ysc(r(0 0.6))
graph export ./images/length_stints_post_y1_data.pdf, replace as(pdf)

reshape wide
***********************
* Length Stints Delta *
***********************
reshape long length_stints_delta, i(years) j(data_type)
label define length_stints_delta_lab 1 "One-Time Entrants" 2 "No Permanent One-Time Entrants" 3 "Re-Entrants First" 4 "Re-Entries"
label values data_type length_stints_delta_lab
graph bar length_stints_delta, over(data_type) over(years) asyvars ytitle("Percent") b1title("Years") legend(position(0) bplacement(nwest)) bar(1, color(51 34 136)) bar(2, color(17 119 51)) bar(3, color(204 102 119)) bar(4, color(82 137 199)) ysc(r(0 0.6))
graph export ./images/length_stints_delta.pdf, replace as(pdf)

reshape wide
**************************
* Length Stints No Delta *
**************************
reshape long length_stints_nodelta, i(years) j(data_type)
label define length_stints_nodelta_lab 1 "One-Time Entrants" 2 "No Permanent One-Time Entrants" 3 "Re-Entrants First" 4 "Re-Entries"
label values data_type length_stints_nodelta_lab
graph bar length_stints_nodelta, over(data_type) over(years) asyvars ytitle("Percent") b1title("Years") legend(position(0) bplacement(nwest)) bar(1, color(51 34 136)) bar(2, color(17 119 51)) bar(3, color(204 102 119)) bar(4, color(82 137 199)) ysc(r(0 0.6)) ytick(0 0.1 0.2 0.3 0.4 0.5 0.6) ylab(0 0.1 0.2 0.3 0.4 0.5 0.6)
graph export ./images/length_stints_nodelta.pdf, replace as(pdf)

reshape wide
rename length_stints_data4 reentry_stints1
rename length_stints_delta4 reentry_stints2
rename length_stints_nodelta4 reentry_stints3
reshape long reentry_stints, i(years) j(data_type)
label define length_stints_lab 1 "Data" 2 "Export Capital" 3 "Sunk Cost"
label values data_type length_stints_lab
graph bar reentry_stints if years < 11, over(data_type) over(years) asyvars ytitle("Percent") b1title("Years") legend(position(0) bplacement(neast)) bar(1, color(82 137 199)) bar(2, color(136 204 238)) bar(3, color(221 204 119))
graph export ./images/length_reentry_stints_models.pdf, replace as(pdf)

reshape wide
rename length_stints_data3 reentry_first_stints1
rename length_stints_delta3 reentry_first_stints2
rename length_stints_nodelta3 reentry_first_stints3
reshape long reentry_first_stints, i(years) j(data_type)
label define length_first_stints_lab 1 "Data" 2 "Export Capital" 3 "Sunk Cost"
label values data_type length_first_stints_lab
graph bar reentry_first_stints if years < 11, over(data_type) over(years) asyvars ytitle("Percent") b1title("Years") legend(position(0) bplacement(neast)) bar(1, color(204 102 119)) bar(2, color(136 204 238)) bar(3, color(221 204 119))
graph export ./images/length_first_stints_models.pdf, replace as(pdf)

reshape wide
rename length_stints_data1 entry_stints1
rename length_stints_delta1 entry_stints2
rename length_stints_nodelta1 entry_stints3
reshape long entry_stints, i(years) j(data_type)
label define length_stints_lab 1 "Data" 2 "Export Capital" 3 "Sunk Cost"
label values data_type length_first_stints_lab
graph bar entry_stints, over(data_type) over(years) asyvars ytitle("Percent") b1title("Years") legend(position(0) bplacement(nwest)) bar(1, color(51 34 136)) bar(2, color(136 204 238)) bar(3, color(221 204 119))
graph export ./images/length_entry_stints_models.pdf, replace as(pdf)

reshape wide
rename length_stints_data2 entry_noperm_stints1
rename length_stints_delta2 entry_noperm_stints2
rename length_stints_nodelta2 entry_noperm_stints3
reshape long entry_noperm_stints, i(years) j(data_type)
label define length_stints_lab 1 "Data" 2 "Export Capital" 3 "Sunk Cost"
label values data_type length_first_stints_lab
graph bar entry_noperm_stints if years < 12, over(data_type) over(years) asyvars ytitle("Percent") b1title("Years") legend(position(0) bplacement(neast))  bar(1, color(17 119 51)) bar(2, color(136 204 238)) bar(3, color(221 204 119))
graph export ./images/length_entry_noperm_stints_models.pdf, replace as(pdf)

**************************************
*           Survival Rates           *
**************************************
clear
svmat survival_delta

label variable survival_delta1 "One-Time Entrants"
label variable survival_delta2 "No Permanent One-Time Entrants"
label variable survival_delta3 "Re-Entrants First Entry"
label variable survival_delta4 "Re-Entries"

svmat survival_nodelta

label variable survival_nodelta1 "One-Time Entrants"
label variable survival_nodelta2 "No Permanent One-Time Entrants"
label variable survival_nodelta3 "Re-Entrants First Entry"
label variable survival_nodelta4 "Re-Entries"

matrix define survival_data = (0.76, 0.61, 0.6, 0.57 \ 0.88, 0.75, 0.7, 0.63 \ 0.92, 0.80, 0.73, 0.87 \ 0.93, 0.80, 0.7, 0.92 \ 0.92, 0.75, 0.54, 0.89)
svmat survival_data

label variable survival_data1 "One-Time Entrants"
label variable survival_data2 "No Permanent One-Time Entrants"
label variable survival_data3 "Re-Entrants First Entry"
label variable survival_data4 "Re-Entries"

gen years_exporting = _n
***********************
* Survival Rates Data *
***********************
label variable years_exporting "Years Exporting"
graph twoway line (survival_data1 survival_data2 survival_data3 survival_data4) years_exporting, ytitle("Survival Rate") legend(position(0) bplacement(south)) lwidth(thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "82 137 199")
graph export ./images/survival_rate_data.pdf, replace as(pdf)

************************
* Survival Rates Delta *
************************
label variable years_exporting "Years Exporting"
graph twoway line (survival_delta1 survival_delta2 survival_delta3 survival_delta4) years_exporting, ytitle("Survival Rate") legend(position(0) bplacement(swest)) lwidth(thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "82 137 199")
graph export ./images/survival_rate_delta.pdf, replace as(pdf)

***************************
* Survival Rates No Delta *
***************************
label variable years_exporting "Years Exporting"
graph twoway line (survival_delta1 survival_delta2 survival_nodelta3 survival_nodelta4) years_exporting, ytitle("Survival Rate") legend(position(0) bplacement(swest)) lwidth(thick thick thick thick) lcolor("51 34 136" "17 119 51" "204 102 119" "82 137 199")
graph export ./images/survival_rate_nodelta.pdf, replace as(pdf)

label variable survival_delta1 "Export Capital"
label variable survival_delta2 "Export Capital"
label variable survival_delta3 "Export Capital"
label variable survival_delta4 "Export Capital"

label variable survival_nodelta1 "Sunk Cost"
label variable survival_nodelta2 "Sunk Cost"
label variable survival_nodelta3 "Sunk Cost"
label variable survival_nodelta4 "Sunk Cost"

label variable survival_data1 "Data"
label variable survival_data2 "Data"
label variable survival_data3 "Data"
label variable survival_data4 "Data"
graph twoway line (survival_delta1 survival_nodelta1 survival_data1) years_exporting, ytitle("Survival Rate") legend(position(0) bplacement(seast)) lwidth(thick thick thick) lcolor("136 204 238" "221 204 119" "51 34 136")
graph export ./images/survival_rate_entrants_models.pdf, replace as(pdf)
graph twoway line (survival_delta2 survival_nodelta2 survival_data2) years_exporting, ytitle("Survival Rate") legend(position(0) bplacement(seast)) lwidth(thick thick thick) lcolor("136 204 238" "221 204 119" "17 119 51")
graph export ./images/survival_rate_entrants_noperm_models.pdf, replace as(pdf)
graph twoway line (survival_delta3 survival_nodelta3 survival_data3) years_exporting, ytitle("Survival Rate") legend(position(0) bplacement(seast)) lwidth(thick thick thick) lcolor("136 204 238" "221 204 119" "204 102 119")
graph export ./images/survival_rate_reentrants_first_models.pdf, replace as(pdf)
graph twoway line (survival_delta4 survival_nodelta4 survival_data4) years_exporting, ytitle("Survival Rate") legend(position(0) bplacement(seast)) lwidth(thick thick thick) lcolor("136 204 238" "221 204 119" "82 137 199")
graph export ./images/survival_rate_reentrants_models.pdf, replace as(pdf)