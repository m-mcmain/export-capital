*! version 1.4.0        9 Dec 2004        Mark Stewart
program define redprob
	version 8
	syntax varlist [if] [in] [ , I(string) T(string) Quadrat(integer 6) FROM(string) ]
	gettoken y 0 : 0
	gettoken xvars 0 : 0 , parse("(")
	gettoken zspec 0 : 0 , parse(",")
	gettoken zvars : zspec , match(parns)
	di _n "Dependent variable = " `"`y'"' _n "x-variables = " `"`xvars'"' _n /*
	*/ "z-variables = " `"`zvars'"'
	xt_iis `i'
	local i "$S_1"
	xt_tis `t'
	local t "$S_1"
	tempvar touse x w
	mark `touse' `if' `in'
	markout `touse' `y' `i' `t'
/* NB: no markout on xvars and zvars */
	local xv1 = word("`xvars'",1)
	di "Assumed lagged dependent variable = " `"`xv1'"'

/* Check to see if outcome varies. */
	quietly count if `touse'
	local n = _result(1)
	quietly count if `y'==0 & `touse'
	local n0 = _result(1)
	if `n0'==0 | `n0'==`n' {
		di _n in blu "outcome does not vary"
		exit
	}

/* Sort data. */
	sort `touse' `i' `t'

/* Get points and weights for Gaussian-Hermite quadrature. */
	ghquad double(`x' `w'), n(`quadrat')

/* Set up macros for ml function. */
	global S_sample "`touse'"
	global S_ivar   "`i'"
	global S_tvar   "`t'"
	global S_x      "`x'"
	global S_w      "`w'"
	global S_quad   "`quadrat'"

/* Set up initial values. */
	tempname llrho0 b0 b00 b1 lllast b ll V b2
	if "`from'"=="" {
	di _n in gr "Pooled Probit Model for t>1"
	probit `y' `xvars' if `touse' & `t'>1
	scalar `llrho0' = _result(2)
	matrix `b0' = get(_b)
	matrix coleq `b0' = `y'

	di _n in gr "Probit Model for t=1"
	probit `y' `zvars' if `touse' & `t'==1
	scalar `llrho0' = `llrho0'+_result(2)
	matrix `b00' = get(_b)
	matrix coleq `b00' = rfper1
	matrix `b0' = `b0' , `b00'

	matrix `b1' = (-0.5, 0)
	matrix colnames `b1' = logitrho:_cons ltheta:_cons
	matrix `b0' = `b0' , `b1'
	}
	else {
	mat `b0' = `from'
	}

/* Set up ml commands. */
	di _n in gr "Iterations for full ML estimation"
	ml model d0 redpmod_ll (`y': `y' = `xvars') (rfper1: `zvars') (logitrho:) (ltheta:) if `touse', /*
	*/ miss nopreserve title("Random-Effects Dynamic Probit Model") init(`b0') search(off) maximize trace grad

	ml display, neq(2) plus
	_diparm logitrho, prob
	_diparm ltheta, prob
	di in smcl in gr "{hline 13}{c +}{hline 64}"
	_diparm logitrho, label("rho") ilogit prob
	_diparm ltheta, label("theta") exp prob
	di in smcl in gr "{hline 13}{c BT}{hline 64}"
	
/* Compute LR test for rho = 0. */
	if "`from'"=="" {
	local chi = 2*(e(ll) - `llrho0')
	#delimit ;
	di in gr "LR test of rho = 0:   chi2(" in ye "1" in gr ")     = "
	   in ye %7.2f `chi' _n
	   in gr "                      Prob > chi2 = " in ye %7.4f
	   chiprob(1,`chi') ;
	#delimit cr
	}

end

/*
	Programs for weights and points for Gaussian-Hermite quadrature follow
	These are due to Bill Sribney at Stata Corp.
*/

* version 1.0.1  29jun1995
program define ghquad
	version 4.0
	local varlist "req new min(2) max(2)"
	local options "N(integer 10)"
	parse "`*'"
	parse "`varlist'", parse(" ")
	local x "`1'"
	local w "`2'"
	if `n' + 2 > _N  {
		di in red  /*
		*/ "`n' + 2 observations needed to compute quadrature points"
		exit 2001
	}
	tempname xx ww
	local i 1
	local m = int((`n' + 1)/2)
	while `i' <= `m' {
		if `i' == 1 {
			scalar `xx' = sqrt(2*`n'+1)-1.85575*(2*`n'+1)^(-1/6)
		}
		else if `i' == 2 { scalar `xx' = `xx'-1.14*`n'^0.426/`xx' }
		else if `i' == 3 { scalar `xx' = 1.86*`xx'-0.86*`x'[1] }
		else if `i' == 4 { scalar `xx' = 1.91*`xx'-0.91*`x'[2] }
		else { scalar `xx' = 2*`xx'-`x'[`i'-2] }
		hermite `n' `xx' `ww'
		qui replace `x' = `xx' in `i'
		qui replace `w' = `ww' in `i'
		local i = `i' + 1
	}
	if mod(`n', 2) == 1 { qui replace `x' = 0 in `m' }
	qui replace `x' = -`x'[`n'+1-_n] in `i'/`n'
	qui replace `w' =  `w'[`n'+1-_n] in `i'/`n'
end

program define hermite  /* integer n, scalar x, scalar w */
	version 4.0
	local n "`1'"
	local x "`2'"
	local w "`3'"
	local last = `n' + 2
	tempvar p
	tempname i
	qui gen double `p' = .
	scalar `i' = 1
	while `i' <= 10 {
		qui replace `p' = 0 in 1
		qui replace `p' = _pi^(-0.25) in 2
		qui replace `p' = `x'*sqrt(2/(_n-2))*`p'[_n-1] /*
		*/	- sqrt((_n-3)/(_n-2))*`p'[_n-2] in 3/`last'
		scalar `w' = sqrt(2*`n')*`p'[`last'-1]
		scalar `x' = `x' - `p'[`last']/`w'
		if abs(`p'[`last']/`w') < 3e-14 {
			scalar `w' = 2/(`w'*`w')
			exit
		}
		scalar `i' = `i' + 1
	}
	di in red "hermite did not converge"
	exit 499
end
	