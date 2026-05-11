*! version 1.4.0        9 Dec 2004        Mark Stewart
program define redpmod_ll /* b log_likelihood */
	version 8
	args todo b f
	local y    = trim("$ML_y1")
	local doit "$S_sample"
	local i    "$S_ivar"
	local t    "$S_tvar"
	tempname rho s2su beta piv u theta u2
	tempvar xb F p zp
	local col = colnumb(`b',"logitrho:_cons")
	scalar `s2su' = sqrt(2*exp(`b'[1,`col']))
	local col = `col'+1
	scalar `theta' = exp(`b'[1,`col'])
	matrix `beta' = `b'[1,"`y':"]
	local k1 = colsof(`beta')+1
	local k2 = colsof(`b')-2
	matrix `piv' = `b'[1,`k1'..`k2']

	quietly {
		matrix score double `xb' = `beta' if `doit' & `t'~=1
		matrix score double `zp' = `piv' if `doit' & `t'==1
		gen double `F' = . in 1
		by `doit' `i': gen double `p' = cond(_n==_N,0,.) if `doit'

		forvalues m = 1/$S_quad {
			scalar `u' = `s2su'*$S_x[`m']
			scalar `u2' = `theta'*`u'
			#delimit ;
			by `doit' `i': replace `F' =
		    	    cond(_n==1,
				cond(`y', normprob(`zp' + `u2'),
	    		      	      1 - normprob(`zp' + `u2')),
				cond(`y', normprob(`xb' + `u'),
	    		      	      1 - normprob(`xb' + `u'))*`F'[_n-1])
		    	    if `doit' ;
			#delimit cr
			replace `p' = `p' + $S_w[`m']*`F' if `doit'
		}

		replace `F' = sum(log(`p'/sqrt(_pi))) if `doit'
		scalar  `f' = `F'[_N]
	}
end
