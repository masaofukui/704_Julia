cd "/Users/fukui/Dropbox (Personal)/teaching/704_2025/704_Julia/Topic2"
capture mkdir ./figure/
clear
import fred OPHNFB,clear

gen yq = qofd(daten)
format yq %tq
tsset yq

keep if yq <= yq(2019,4)
gen log_lp = log(OPHNFB)
gen dllp = log_lp - L.log_lp
tw line dllp yq
sum dllp
gen d_dllp = dllp - `r(mean)'

tsfilter bk llp_bk = log_lp

reg llp_bk L.llp_bk
global rho_z = (_b[L1.llp_bk])^(1/4)
di $rho_z
tw scatter llp_bk L.llp_bk


gen ym = mofd(daten)
format ym %tm
tsset ym
tsfill
di $rho_z

gen llp_bk_int = $rho_z*L.llp_bk
replace llp_bk_int = llp_bk if llp_bk_int ==.
replace llp_bk_int = $rho_z*L.llp_bk_int if llp_bk_int ==.
gen llp_bk_shock = llp_bk_int - $rho_z*L.llp_bk_int
replace llp_bk_shock = llp_bk_int if llp_bk_shock == .

drop if llp_bk_shock == .
outsheet daten ym llp_bk_shock using ./data/labor_productivity_bk.csv,replace comma
save ./data/bk_prod,replace

/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 
* simulate in Julia
/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 
insheet using "./data/u_sim.csv",clear
rename v1 simualted_u
gen nn = _n
save ./data/simualted_u,replace


/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 
import fred UNRATE, clear
gen ym = mofd(daten)
save ./data/unrate_save,replace
/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 
use ./data/bk_prod,clear
gen nn = _n
merge 1:1 nn using ./data/simualted_u,nogen keep(1 3)
merge m:1 ym using ./data/unrate_save,nogen keep(1 3)
replace simualted_u = simualted_u*100


*nbercycles llp_bk, file("bk_filter.do")
local yms = ym(1950,1)
local ymend = ym(2019,10)
twoway function y=.0288060759007931,range(-78 -68) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(-29 -21) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(3 13) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(119 130) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(166 182) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(240 246) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(258 274) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(366 374) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(494 502) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(575 593) recast(area) color(gs14) base(-.0242569060809911) || /// 
tsline llp_bk , xlabel(`yms'(120)`ymend',format(%tmCCYY)) legend(off) graphregion(color(white)) ///
xsize(10) ytitle("Detrended log Labor Productivity") scheme(s2color)
graph export ./figure/labor_productivity.pdf,replace

local yms = ym(1950,1)
local ymend = ym(2019,10)
twoway function y= 12,range(-78 -68) recast(area) color(gs14) base(2)|| /// 
function y=12,range(-29 -21) recast(area) color(gs14) base(2) || /// 
function y=12,range(3 13) recast(area) color(gs14) base(2) || /// 
function y=12,range(119 130) recast(area) color(gs14) base(2) || /// 
function y=12,range(166 182) recast(area) color(gs14) base(2) || /// 
function y=12,range(240 246) recast(area) color(gs14) base(2) || /// 
function y=12,range(258 274) recast(area) color(gs14) base(2) || /// 
function y=12,range(366 374) recast(area) color(gs14) base(2) || /// 
function y=12,range(494 502) recast(area) color(gs14) base(2) || /// 
function y=12,range(575 593) recast(area) color(gs14) base(2) || /// 
tsline simualted_u , lw(0.7) lc(maroon) || ///
( tsline UNRATE , lw(0.7) lc(navy) lpattern(dash) ) , xlabel(`yms'(120)`ymend',format(%tmCCYY)) graphregion(color(white)) ///
xsize(10) ytitle("Unemployment Rate (%)") ///
legend(order(11 "Simulated u" 12 "Actual u")) scheme(s2color)
graph export ./figure/shimer_puzzle.pdf,replace

sum simualted_u
sum UNRATE


tsfilter bk urate_bk = UNRATE
tsline urate_bk

local yms = ym(1950,1)
local ymend = ym(2019,10)
twoway function y=.0288060759007931,range(-78 -68) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(-29 -21) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(3 13) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(119 130) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(166 182) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(240 246) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(258 274) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(366 374) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(494 502) recast(area) color(gs14) base(-.0242569060809911) || /// 
function y=.0288060759007931,range(575 593) recast(area) color(gs14) base(-.0242569060809911) || /// 
(tsline llp_bk, lw(0.7)) (tsline UNRATE, lw(0.7) yaxis(2) lc(maroon) lp(dash)) , xlabel(`yms'(120)`ymend',format(%tmCCYY))  graphregion(color(white)) ///
xsize(10) ytitle("Detrended log Labor Productivity") ytitle("Detrended Unemployment rate",axis(2)) legend(order(11 "Labor Productivity" 12 "Unemployment Rate")) ///
xtitle("") scheme(s2color)
graph export ./figure/labor_productivity_urate.pdf,replace



/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 
import excel using "./original_data/ie_data.xls", sheet("Data") clear cellrange(A8) first
keep A H
gen year = floor(A)
gen month = mod(A,1)*100
gen ym = ym(year,month)
rename H real_P
keep ym real_P
format ym %tm

merge m:1 ym using ./data/unrate_save, nogen keep(1 3)
save ./data/realP,replace

use ./data/realP,clear
tsset ym
gen lp = log(real_P)
tsfilter bk urate_bk = UNRATE
tsfilter bk lp_bk = lp

reg lp ym
predict lp_l, residual

local yms = ym(1950,1)
local ymend = ym(2019,10)
keep if inrange(ym,`yms',`ymend')
twoway function y=1,range(-78 -68) recast(area) color(gs14) base(-1) || /// 
function y=1,range(-29 -21) recast(area) color(gs14)  base(-1) || /// 
function y=1,range(3 13) recast(area) color(gs14)  base(-1) || /// 
function y=1,range(119 130) recast(area) color(gs14)  base(-1) || /// 
function y=1,range(166 182) recast(area) color(gs14)  base(-1) || /// 
function y=1,range(240 246) recast(area) color(gs14)  base(-1) || /// 
function y=1,range(258 274) recast(area) color(gs14)  base(-1) || /// 
function y=1,range(366 374) recast(area) color(gs14)  base(-1) || /// 
function y=1,range(494 502) recast(area) color(gs14)  base(-1) || /// 
function y=1,range(575 593) recast(area) color(gs14)  base(-1) || /// 
(tsline lp_l, lw(0.7)) (tsline UNRATE, lw(0.7) yaxis(2) lc(maroon) lp(dash)) , xlabel(`yms'(120)`ymend',format(%tmCCYY))  graphregion(color(white)) ///
xsize(10) ytitle("Detrended log Real Stock Price (S&P 500)") ytitle("Detrended Unemployment rate",axis(2)) legend(order(11 "Stock Price" 12 "Unemployment Rate")) ///
xtitle("") scheme(s2color)
graph export ./figure/stock_price_urate.pdf,replace










