cd "/Users/fukui/Dropbox (Personal)/teaching/704_2025/704_Julia/Topic6/data/"
wbopendata, language(en - English) indicator(FS.AST.PRVT.GD.ZS; NY.GDP.MKTP.CN) long clear
keep countrycode year fs_ast ny_gdp
save ./working_data/credit_gdp,replace



/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 
use ./original_data/pwt1001.dta,clear
merge m:1 countrycode year using ./working_data/credit_gdp, nogen keep(1 3)
gen K_Y = rgdpna/rnna
gen ltfp = log(rtfpna)
gen lcwtfp = log(cwtfp)
gen lctfp = log(ctfp)

gen D_K = K_Y*fs_ast/100

tw scatter lctfp fs_ast if year == 2019, mlabel(countrycode)
tw ( scatter lcwtfp fs_ast if year == 2019 & fs_ast!=., mlabel(countrycode) mlabposition(0) msymbol(none) ) ///
	( lfit lcwtfp fs_ast if year == 2019 & fs_ast!=., lw(0.8)), legend(off) name(welfare,replace) ///
	xtitle("Credit to GDP Ratio in 2019 (%)") ytitle("log TFP in 2019 (USA=0)") xlab(,nogrid) 
graph export ./figure/credit_gdp_tfp.pdf

tw ( scatter cwtfp D_K if year == 2019 & D_K!=., mlabel(countrycode) mlabposition(0) msymbol(none) ) ///
	( lfit cwtfp D_K if year == 2019 & D_K!=., lw(0.8)) , ///
	 legend(off) name(welfare,replace) ///
	xtitle("Credit to GDP Ratio in 2019 (%)") ytitle("log TFP in 2019 (USA=0)") xlab(,nogrid) ///


tw ( scatter cwtfp fs_ast if year == 2019 & fs_ast!=., mlabel(countrycode) ) ///
	( lfit cwtfp fs_ast if year == 2019 & fs_ast!=.), legend(position(6)) name(tfp,replace)

tw line fs_ast year if countrycode == "USA"


/*  ------- break line ------- (by the Stata editor for macOS (piu_sign) )  */ 

use "./original_data/MV_credit_data.dta",clear
rename iso3c countrycode
merge m:1 countrycode year using ./working_data/credit_gdp, nogen keep(1 3)
merge m:1 countrycode year using ./original_data/pwt1001.dta, nogen keep(1 3)
replace ny_gdp = ny_gdp/1000000
gen total_gdp = Total/ny_gdp
gen hh_gdp = HH/ny_gdp
gen firm_gdp = Corp/ny_gdp
tw (line total_gdp year ) (line hh_gdp year ) (line firm_gdp year ) if countrycode == "USA"
tw line hh_gdp year if countrycode == "USA"

tw ( scatter cwtfp total_gdp if year == 2014 & total_gdp!=., mlabel(countrycode) ) ///
	( lfit cwtfp total_gdp if year == 2014 & total_gdp!=.) ///
	if total_gdp <= 3, legend(position(6)) name(tfp,replace)

tw ( scatter lcwtfp firm_gdp if year == 2014 & total_gdp!=., mlabel(countrycode) ) ///
	( lfit lcwtfp firm_gdp if year == 2014 & total_gdp!=.) ///
	if firm_gdp <= 3, legend(position(6)) name(tfp,replace)


tw ( scatter cwtfp hh_gdp if year == 2014 & total_gdp!=., mlabel(countrycode) ) ///
	( lfit cwtfp hh_gdp if year == 2014 & total_gdp!=.) ///
	if hh_gdp <= 3, legend(position(6)) name(tfp,replace)


use "/Users/fukui/Dropbox (Personal)/FTY_Forward_Guidance_Puzzle/data/Jorda_dataset/JSTdatasetR6.dta",clear
rename iso countrycode 
gen tloan_gdp = tloan/gdp
gen thh_gdp = thh/gdp
gen tbus_gdp = tbus/gdp

merge m:1 countrycode year using ./original_data/pwt1001.dta, nogen keep(1 3)

tw line tloan_gdp year if countrycode== "JPN"

tw (line thh_gdp year, lw(0.8)) (line tbus_gdp year, lw(0.8))  if countrycode== "USA" & !missing(thh_gdp) ///
, legend(order(1 "Household debt" 2 "Firm debt") rows(1) position(6)) ytitle("% of GDP") xtitle("Year") 
graph export ./figure/hh_debt.pdf

tw ( scatter cwtfp thh_gdp if year == 2014 & thh_gdp!=., mlabel(countrycode) ) ///
	( lfit cwtfp thh_gdp if year == 2014 & thh_gdp!=.) ///
, legend(position(6)) name(tfp,replace)

