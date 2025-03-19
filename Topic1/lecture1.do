cd "/Users/fukui/Dropbox (Personal)/teaching/704_2024/stata"
capture mkdir ./figure/
clear
import fred CIVPART LNS11300002 LNS11300001 LNS11300060 LRAC25MAUSM156S LRAC25FEUSM156S
global lwset = 0.8

forval yr = 1950(10)2020{
  local yy = mdy(1,1,`yr')
  local yrlabel `yrlabel' `yy'
}

tw ( line LNS11300001  daten , lw(${lwset}) ) ///
 ( line CIVPART daten , lw(${lwset}) lc(forest_green) ) ///
  ( line LNS11300002 daten ,  lw(${lwset})lc(maroon) ),  ///
  xlabel(`yrlabel',format(%tdCCYY)) ///
  graphregion(color(white)) xtitle("") ytitle("%") ///
  legend(order(2 "Total" 1 "Male" 3 "Female") rows(1) position(6)) ///
   title("Labor Force Participation Rate") xsize(9) ysize(5) xlab(,nogrid)
graph export ./figure/lfp.pdf,replace



gen ym = mofd(daten)
format ym %tm
tsset ym

*nbercycles CIVPART if tin(1980m1,2023m2) , file(lfp.do) replace


forval yr = 1980(10)2020{
  local yy = ym(`yr',1)
  local yrlabel `yrlabel' `yy'
}
local yms = ym(1980,1)
twoway function y=67.9730030822754,range(240 246) recast(area) color(gs14) base(59.49899848937987) || /// 
function y=67.9730030822754,range(258 274) recast(area) color(gs14) base(59.49899848937987) || /// 
function y=67.9730030822754,range(366 374) recast(area) color(gs14) base(59.49899848937987) || /// 
function y=67.9730030822754,range(494 502) recast(area) color(gs14) base(59.49899848937987) || /// 
function y=67.9730030822754,range(575 593) recast(area) color(gs14) base(59.49899848937987) || /// 
function y=67.9730030822754,range(721 723) recast(area) color(gs14) base(59.49899848937987) || /// 
(line CIVPART ym , lw(${lwset}) lc(navy)) if ym >=`yms' ,  ///
 legend(order(7 1 "Recession")) graphregion(color(white)) legend(off) ytitle("%") title("Labor Force Participation Rate") ///
 xsize(10) ysize(6) xlabel(`yrlabel',format(%tmCCYY))
 graph export ./figure/lfp_total.pdf,replace


forval yr = 1950(10)2020{
  local yy = mdy(1,1,`yr')
  local yrlabel `yrlabel' `yy'
}

tw ( line LRAC25MAUSM156S  daten , lw(${lwset}) ) ///
 ( line LNS11300060 daten , lw(${lwset}) lc(forest_green) ) ///
  ( line LRAC25FEUSM156S daten ,  lw(${lwset})lc(maroon) ),  ///
  xlabel(`yrlabel',format(%tdCCYY)) ///
  graphregion(color(white)) xtitle("") ytitle("%") ///
  legend(order(2 "Total" 1 "Male" 3 "Female") rows(1)) title("Labor Force Participation Rate: 25-54 yrs old") xsize(9) ysize(5) ///
  ylabel(20(20)100) ///
  legend(position(6)) xlab(,nogrid)
graph export ./figure/lfp_prime_age.pdf,replace



/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 
clear
global lwset = 0.8

import fred LNS17400000 LNS17100000 LNS17200000 LNS17900000  LNS17800000 LNS17600000 USRECM,clear
keep if daten >= mdy(1,1,1990)

sum LNS17400000
foreach var of varlist LN*{
  replace `var' = `var'/1000
}

global covid = 1
label var LNS17400000 "E to U"
rename LNS17400000 etou_f
label var LNS17100000 "U to E"
rename LNS17100000 utoe_f
label var LNS17200000 "N to E"
rename LNS17200000 ntoe_f
label var LNS17900000 "U to N"
rename LNS17900000 uton_f
label var LNS17800000 "E to N"
rename LNS17800000 eton_f
label var LNS17600000 "N to U"
rename LNS17600000 ntou_f


foreach var of varlist *_f {
  local vtext : variable label `var' 

  if inlist("`var'", "ntoe_f", "eton_f"){
    local recnum = 6
  }
  else{
        local recnum = 4
  }
  if $covid == 1{
    replace USRECM = 20*(USRECM > 0)
    local timecon
  }
  else{
    replace USRECM = `recnum'*(USRECM > 0)
    local timecon "if daten <= mdy(12,31,2019)"
  }

  tw (area USRECM daten,color(gs14))  ///
  (line `var' daten, lw(${lwset}) color(maroon)) ///
  `timecon', ///
  ytitle("Millions of Persons") ///
  graphregion(color(white)) title("`vtext'") ///
   xlabel(,format(%tdCCYY))  xtitle("") legend(off) name(`var',replace) 
}
graph combine etou_f utoe_f ntoe_f eton_f uton_f ntou_f, ///
graphregion(color(white)) xsize(10) ysize(5) rows(3)
graph export ./figure/Flow_all_covid_$covid.pdf,replace


drop datestr USRECM
save ./data/fred_labor_flow,replace


/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 

infix year 1-4 month 8-9 unrate 11-20 using "./original_data/m08292a.dat",clear
save ./data/urate_1929_40,replace
infix year 1-4 month 8-9 unrate 11-20 using "./original_data/m08292b.dat",clear
save ./data/urate_1940_46,replace
infix year 1-4 month 8-9 unrate 11-20 using "./original_data/m08292c.dat",clear
save ./data/urate_1947_66,replace

use ./data/urate_1929_40,clear
drop if year>=1940
append using ./data/urate_1940_46
append using ./data/urate_1947_66
gen ym = ym(year,month)
drop if year >= 1948
keep ym unrate
save ./data/unrate_old_all,replace


import fred unrate USRECM,clear
gen ym = mofd(daten)
rename UNRATE unrate
keep unrate ym USRECM
format ym %tm
append using ./data/unrate_old_all
save ./data/unrate_all,replace
keep if ym >= ym(1929,1)
sort ym
replace USRECM = 26*USRECM


forval yr = 1930(10)2020{
  local yy = ym(`yr',1)
  local yrlabel `yrlabel' `yy'
}
tw (area USRECM ym,color(gs14))  ///
(line unrate ym, lc(navy) lw(0.7)), ///
yscale(range(0,26)) graphregion(color(white)) ///
xlabel(`yrlabel',format(%tmCCYY)) xtitle("") ///
legend(order(1 "NBER recession" 2 "Unemployment rate") position(6) rows(1)) ///
title("Unemployment Rate") ytitle("%") ///
xsize(12) ysize(6) note("Data: NBER Macro History Database and CPS")
graph export ./figure/unrate.pdf,replace

import fred LNS14000001 LNS14000002,clear



forval yr = 1950(10)2020{
  local yy = mdy(1,1,`yr')
  local yrlabel `yrlabel' `yy'
}
tw  ( line LNS14000001 daten , lw(${lwset}) lc(navy) ) ///
  ( line LNS14000002 daten ,  lw(${lwset})lc(maroon) lp(dash) ),  ///
  xlabel(`yrlabel',format(%tdCCYY)) ///
  graphregion(color(white)) xtitle("") ytitle("%") ///
  legend(order(1 "Male" 2 "Female") rows(1)) title("Unemployment Rate by Gender") xsize(9) ysize(5) ///
  legend(position(6)) xlab(,nogrid)
graph export ./figure/unrate_by_gender.pdf,replace

sum LNS14000001
sum LNS14000002


/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 

use "./original_data/cps_00048.dta",clear
drop if cpsidp == 0
gen employed = 0
replace employed = 1 if inlist(empstat,10,12)

gen in_lf = 0
replace in_lf = 1 if labforce == 2
gen ym = ym(year,month)
format ym %tm
format cpsidp %15.0f
xtset cpsidp ym, monthly

gen UtoE = 1 if L.employed == 0 & L.in_lf == 1 & employed == 1
gen EtoU = 1 if L.employed == 1 & in_lf == 1 & employed == 0
gen EtoN = 1 if L.employed == 1 & in_lf == 0 & employed == 0
gen UtoN = 1 if L.employed == 0 & L.in_lf == 1 & in_lf == 0
gen NtoE = 1 if employed == 1 & L.in_lf == 0 
gen NtoU = 1 if employed == 0 & in_lf == 1 & L.in_lf == 0 

collapse (sum) UtoE EtoU EtoN UtoN NtoE NtoU [iw = wtfinl],by(ym)
tw (line UtoE ym) (line EtoU ym)

outsheet using ../R/CPS_aggregate.csv, replace comma


tset ym
sort ym
foreach var of varlist *{
  replace `var' =  L.`var' if `var' <= 206858 & "`var'" != "ym"
}
outsheet using ../R/CPS_aggregate_nonzero.csv, replace comma


/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 

import fred USRECM,clear
save ./data/recession,replace

insheet using "../R/adjusted.csv",clear
global covid = 1
gen daten =date(ym,"MY")
format daten %td

merge m:1 daten using ./data/recession, nogen keep(1 3)
merge m:1 daten using ./data/fred_labor_flow, nogen keep(1 3)

foreach var of varlist utoe etou ntou ntoe eton uton{
  replace `var' = `var'/1000000
  gen r_`var'_f = `var'_f/`var'
  sum r_`var'_f if abs(r_`var'_f) < 2
  gen `var'_fc = `var'_f
  replace `var'_fc = `var'*`r(mean)' if `var'_fc == .
  gen `var'_check = `var'*`r(mean)'

}
tw (line ntou_check daten) (line ntou daten ) (line ntou_fc daten )
label var utoe_fc "U to E"
label var etou_fc "E to U"
label var ntou_fc "N to U"
label var ntoe_fc "N to E"
label var eton_fc "E to N"
label var uton_fc "U to N"


foreach var of varlist *_fc {
  local vtext : variable label `var' 
  sum `var'
  replace `var' = `var'[_n-1] if `var' <= `r(mean)'/3

  if ${covid} == 0{
    if inlist("`var'", "ntoe_fc", "eton_fc"){
      local recnum = 6
    }
    else{
      local recnum = 4
    }
    replace USRECM = `recnum'*(USRECM>0)
        local tcond "if daten <= mdy(12,31,2019)"
  }
  else{
    if inlist("`var'", "etou_fc"){
      local recnum = 20
    }
    else if inlist("`var'", "utoe_fc","eton_fc"){
      local recnum = 10
    }
    else if inlist("`var'", "uton_fc","ntou_fc"){
      local recnum = 5
    }
    else if inlist("`var'", "ntoe_fc"){
      local recnum = 6
    }

    replace USRECM = `recnum'*(USRECM>0)
    local tcond ""
  }
  tw (area USRECM daten,color(gs14))  ///
  (line `var' daten, lw(${lwset})) ///
  `tcond', ///
  ytitle("Millions of Persons") ///
  graphregion(color(white)) title("`vtext'") ///
   xlabel(,format(%tdCCYY))  xtitle("") legend(off) name(`var', replace) 
}

graph combine etou_fc utoe_fc ntoe_fc eton_fc uton_fc ntou_fc, graphregion(color(white)) rows(3) ///
xsize(10) ysize(5)
graph export ./figure/Flow_CPS_all_covid_$covid.pdf,replace

gen de = utoe_fc - etou_fc
gen dn = uton_fc - ntou_fc

replace USRECM = `recnum'*(USRECM>0)
gen recup = 0.6*(USRECM>0)
gen reclow = -0.6*(USRECM>0)

drop ym
gen ym = mofd(daten)
format ym %tm
tset ym
do nber2.do
/* nbercycles de if tin(1980m1,2019m12) , file(nber2.do) replace ///
gropt( legend(off) graphregion(color(white)) color(maroon)) 

nbercycles dn if tin(1980m1,2019m12) , file(nber3.do) replace ///
gropt( legend(off) graphregion(color(white)) color(maroon)) 
*/ 
tw (area recup reclow daten,color(gs14)) (line de daten)  if daten <= mdy(12,31,2019)
tw line dn daten if daten <= mdy(12,31,2019)

keep daten *_fc
save ./data/flow,replace

/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 

infix time 1-8 eu_shimer 9-20  using ./original_data/shimer/eu.dat,clear
tw line eu_shimer time 
gen q = (time - floor(time))*4 + 1
gen y = floor(time)
gen yq = yq(y,q)
gen daten = dofq(yq)
format daten %td
keep daten eu_shimer
replace eu_shimer = eu_shimer*100
save ./data/eu_shimer,replace

infix time 1-8 ue_shimer 9-20  using ./original_data/shimer/ue.dat,clear
tw line ue_shimer time 
gen q = (time - floor(time))*4 + 1
gen y = floor(time)
gen yq = yq(y,q)
gen daten = dofq(yq)
format daten %td
keep daten ue_shimer
replace ue_shimer = ue_shimer*100
save ./data/ue_shimer,replace


import delimited  using "./original_data/shimer/find-prob.dat",clear
destring v1, ignore("{") replace
destring v2, ignore("}") replace
rename v2 f_shimer
rename v1 time
tw line f_shimer time 
gen q = (time - floor(time))*4 + 1
gen y = floor(time)
gen yq = yq(y,q)
gen daten = dofq(yq)
format daten %td
keep daten f_shimer
replace f_shimer = f_shimer*100
save ./data/f_shimer,replace


import delimited  using "./original_data/shimer/sep-prob.dat",clear
destring v1, ignore("{") replace
destring v2, ignore("}") replace
rename v2 s_shimer
rename v1 time
tw line s_shimer time 
gen q = (time - floor(time))*4 + 1
gen y = floor(time)
gen yq = yq(y,q)
gen daten = dofq(yq)
format daten %td
keep daten s_shimer
replace s_shimer = s_shimer*100
save ./data/s_shimer,replace



/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */
clear 
import fred CE16OV UNEMPLOY UNRATE
rename CE16OV employed
rename UNEMPLOY unemployed
replace employed = employed/1000
replace unemployed = unemployed/1000
rename UNRATE unrate
merge 1:1 daten using ./data/flow, nogen keep(1 3)
merge m:1 daten using ./data/recession, nogen keep(1 3)

gen ym = mofd(daten)
format ym %tm
tsset ym
gen job_finding_bias = utoe_fc/L.unemployed
gen separation_bias = etou_fc/L.employed

save ./data/fred_unrate,replace
outsheet using ./data/fred_unrate.csv,replace


/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 
* run time_aggreation.R
/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 


insheet using  ./data/fred_unrate_bias_corrected.csv,clear
rename usrecm USRECM
drop daten ym
gen daten = date(datestr,"YMD")
gen ym = mofd(daten)
format daten %td
format ym %tm
save ./data/fred_unrate_bias_corrected,replace
merge m:1 daten using ./data/eu_shimer, nogen keep(1 3)
merge m:1 daten using ./data/ue_shimer, nogen keep(1 3)
merge m:1 daten using ./data/f_shimer, nogen keep(1 3)
merge m:1 daten using ./data/s_shimer, nogen keep(1 3)
keep if inrange(daten,mdy(1,1,1980),mdy(12,31,2019))

sort ym
tsset ym

replace separation = separation*100
replace job_finding = job_finding*100
replace separation_bias = separation_bias*100
replace job_finding_bias = job_finding_bias*100


gen hatu = separation/(separation + job_finding)*100
gen lhatu = log( hatu /(100-hatu) )
gen ls = log(separation)
gen lf = log(job_finding)
gen mlf = -lf

sum mlf
gen lf_demean = mlf - `r(mean)'
sum ls
gen ls_demean = ls - `r(mean)'
sum lhatu
gen lhatu_demean = lhatu - `r(mean)'


local yms = ym(1980,1)
local ymd = ym(2020,1)
tw function y=12,range(240 246) recast(area) color(gs12) base(2) || /// 
function y=12,range(258 274) recast(area) color(gs12) base(2) || /// 
function y=12,range(366 374) recast(area) color(gs12) base(2) || /// 
function y=12,range(494 502) recast(area) color(gs12) base(2) || /// 
function y=12,range(575 593) recast(area) color(gs12) base(2) || /// 
( line hatu ym, lw(0.5) lc(navy)) ( line unrate ym, lw(0.5) lp(dash) lc(maroon)), ///
graphregion(color(white)) legend(order(6 "hat u" 7 "Actual u")) ///
xsize(9) xlabel(`yms'(120)`ymd',format(%tmCCYY)) ytitle("%")
graph export ./figure/hatu_fit.pdf,replace


sum job_finding
gen hatus = separation/(separation + `r(mean)')*100

sum separation
gen hatuf = `r(mean)'/(`r(mean)' + job_finding)*100

tw ( line unrate ym , lw(0.5)) ///
( line hatus ym, lp(dash) lw(0.5) ) ///
( line hatuf ym,lp(dash_dot) lw(0.5) ) ///
, graphregion(color(white)) legend(order(1 "Data" 2 "Only s varying" 3 "Only f varying") rows(1)) ///
xsize(9)  xlabel(,format(%tmCCYY)) xtitle("") ytitle("%")
graph export ./figure/hatu_fit_decom.pdf,replace




tw (function x, lc(forest_green) lp(dash) range(-0.6 0.6)) ///
(scatter lf_demean lhatu_demean,msymbol(oh) color(navy)) ///
(lfit lf_demean lhatu_demean, lc(maroon) lw(0.6)) ///
, xlabel(-0.6(0.2)0.6) graphregion(color(white)) legend(off) ///
ytitle("- log f (demeaned)") xtitle("log (hat u/ ( 1- hat u)) (demeaned)") ///
name(f_sc,replace) title(" - log f vs log u")


tw (function x, lc(forest_green) lp(dash) range(-0.6 0.6)) ///
(scatter lf_demean lhatu_demean,msymbol(oh) color(navy)) ///
(lfit lf_demean lhatu_demean, lc(maroon) lw(0.6)) ///
, xlabel(-0.6(0.2)0.6) graphregion(color(white)) legend(off) ///
ytitle("log s (demeaned)") xtitle("log (hat u/ ( 1- hat u)) (demeaned)") ///
name(s_sc,replace) title(" log s vs log u")

graph combine s_sc f_sc , xsize(10) graphregion(color(white)) iscale(1)
graph export ./figure/decom.pdf,replace

tw ( line separation_bias daten ) ( line separation daten )

replace USRECM = 2.5*(USRECM > 0)
tw (area USRECM daten,color(gs14)) ///
( line separation daten  ) ///
if separation !=. & daten <= mdy(12,31,2019), graphregion(color(white))

do ./separation.do 




reg ls lhatu  if daten <= mdy(1,1,2020)
reg lf lhatu  if daten <= mdy(1,1,2019)
tw (scatter ls lhatu)




sort ym
gen hatu_shimer = eu_shimer/(eu_shimer + ue_shimer)*100
gen lhatu_shimer = log( hatu_shimer /(100-hatu_shimer) )

gen hatu_shimer_fs = s_shimer/(s_shimer + f_shimer)*100
gen lhatu_shimer_fs = log( hatu_shimer_fs /(100-hatu_shimer_fs) )

gen lue_shimer = log(ue_shimer)
gen leu_shimer = log(eu_shimer)
gen lf_shimer = log(f_shimer)
gen ls_shimer = log(s_shimer)

reg lf_shimer lhatu_shimer_fs
reg ls_shimer lhatu_shimer_fs

sum s_shimer
gen hatu_shimer_f = `r(mean)'/(`r(mean)' + f_shimer)*100


tw ( line hatu_shimer_fs daten  ) ( line unrate daten  ) ///
( line hatu_shimer_f daten  ) ( line hatu_shimer_s daten  )
reg S12.lf_shimer S12.lhatu_shimer_fs 

sum ue_shimer
gen hatus_shimer = eu_shimer/(eu_shimer + `r(mean)')*100

sum f_shimer
gen hatu_shimer_s = s_shimer/(s_shimer + `r(mean)')*100

tw ( line hatu_shimer_s daten) ( line hatu_shimer_fs daten)
sum eu_shimer
gen hatuf_shimer = `r(mean)'/(`r(mean)' + ue_shimer)*100

reg unrate hatuf_shimer
reg unrate hatus_shimer

tw (line hatuf_shimer daten) (line hatus_shimer daten) (line unrate daten)

reg S36.lue_shimer S36.lhatu_shimer

reg lhatu_shimer lue_shimer
reg lhatu_shimer leu_shimer
*keep if daten >= mdy(1,1,1980)
sort ym
gen job_finding = utoe_fc/L.unemployed*100
gen separation = etou_fc/L.employed*100

tw (line separation daten )  (line eu_shimer daten ) if daten <= mdy(1,1,2008)
tw (line job_finding daten )  (line ue_shimer daten ) if daten <= mdy(1,1,2008)

replace USRECM = 2.5*(USRECM > 0)
tw (area USRECM daten,color(gs14)) ///
( line separation daten  ) ///
if separation !=. & daten <= mdy(12,31,2019), graphregion(color(white))

do ./separation.do 
/*
nbercycles separation if tin(1980m1,2019m12) , file(separation.do) replace ///
gropt( legend(off) graphregion(color(white)))
*/





/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 
/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */
clear
import fred JTSJOL unrate UNEMPLOY USRECM ,clear
gen ym = mofd(daten)
replace JTSJOL = JTSJOL/1000
save ./fred_v,replace

use ./fred_v,clear
merge m:1 daten using ./data/fred_unrate_bias_corrected, nogen keep(1 3)
replace USRECM = (USRECM > 0)*12
local yms = ym(2000,1)
local ymd = ym(2022,12)
gen lv = log(JTSJOL)
gen lu = log(UNEMPLOY)
replace UNEMPLOY = UNEMPLOY/1000

gen vrate = UNEMPLOY/UNRATE*JTSJOL
gen lf = log(job_finding)
gen lvu = log(vrate/UNRATE)
/*
format ym %tm
tsset ym
nbercycles JTSJOL, file(vacancy.do)
*/

reg lf lvu


forval yr = 2000(5)2025{
  local yy = ym(`yr',1)
  local yrlabel `yrlabel' `yy'
}

tw function y=11.97354953765869,range(494 502) recast(area) color(gs14) base(2.209680111408233) || /// 
function y=11.97354953765869,range(575 593) recast(area) color(gs14) base(2.209680111408233) || /// 
function y=11.97354953765869,range(721 723) recast(area) color(gs14) base(2.209680111408233) || /// 
( line JTSJOL ym, lc(maroon) lw(0.5))  if JTSJOL != .,xlabel(`yrlabel',format(%tmCCYY)) ///
graphregion(color(white)) ytitle("Job Openings (Millions)") xtitle("") legend(off)
graph export ./figure/vacancy.pdf, replace


tw ( scatter vrate UNRATE, ms(oh) ) if ym <= ym(2009,12) ///
, graphregion(color(white)) ytitle("Vacancy rate (%)") xtitle("Unemployment rate (%)") ///
legend(on order(1 "2000-2009")  rows(1)) xlabel(5(5)15) ylabel(5(5)20) ///
 xsize(10) ysize(6) yscale(range(3.5 20)) ///
 legend(position(6)) xlab(,nogrid)
graph export ./figure/uv1.pdf, replace

tw ( scatter vrate UNRATE  if ym <= ym(2009,12), ms(oh) ) ///
 ( scatter vrate UNRATE  if inrange(ym,ym(2010,1),ym(2019,12)), ms(th) ) ///
, graphregion(color(white)) ytitle("Vacancy rate (%)") xtitle("Unemployment rate (%)") ///
legend(order(1 "2000-2009" 2 "2010-2019")  rows(1)) xlabel(5(5)15) ylabel(5(5)20)  ///
 xsize(10) ysize(6) yscale(range(3.5 20)) ///
  legend(position(6)) xlab(,nogrid)
graph export ./figure/uv2.pdf, replace


tw ( scatter vrate UNRATE  if ym <= ym(2009,12), ms(oh) ) ///
 ( scatter vrate UNRATE  if inrange(ym,ym(2010,1),ym(2019,12)), ms(th) ) ///
  ( scatter vrate UNRATE  if inrange(ym,ym(2020,1),ym(2022,12)), ms(sh) ) ///
, graphregion(color(white)) ytitle("Vacancy rate (%)") xtitle("Unemployment rate (%)") ///
legend(order(1 "2000-2009" 2 "2010-2019" 3 "2020-2022") rows(1)) xlabel(5(5)15) ylabel(5(5)20)  ///
 xsize(10) ysize(6) yscale(range(3.5 20)) ///
  legend(position(6)) xlab(,nogrid)
graph export ./figure/uv3.pdf, replace


tw ( scatter vrate UNRATE  if ym <= ym(2009,12), ms(oh) ) ///
 ( scatter vrate UNRATE  if inrange(ym,ym(2010,1),ym(2019,12)), ms(th) ) ///
  ( scatter vrate UNRATE  if inrange(ym,ym(2020,1),ym(2022,12)), ms(sh) ) ///
    ( scatter vrate UNRATE  if ym > ym(2022,12), ms(x) msize(3) ) ///
, graphregion(color(white)) ytitle("Vacancy rate (%)") xtitle("Unemployment rate (%)") ///
legend(order(1 "2000-2009" 2 "2010-2019" 3 "2020-2022" 4 "2022-2025") rows(1)) xlabel(5(5)15) ylabel(5(5)20)  ///
 xsize(10) ysize(6) yscale(range(3.5 20)) ///
  legend(position(6)) xlab(,nogrid)
graph export ./figure/uv4.pdf, replace




tw (scatter lf lvu, ms(oh) ) (lfit lf lvu ,lw(0.5)) ///
if inrange(ym,ym(2000,1),ym(2019,12)) , graphregion(color(white)) ///
xtitle("log (v/u)") ytitle("log f") title("2000-2019") note("Note: The red line is the best linear fit.") ///
legend(off) xsize(7) ///
 legend(position(6)) xlab(,nogrid)

graph export ./figure/lf_lvu.pdf,replace


