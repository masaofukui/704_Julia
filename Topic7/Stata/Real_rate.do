cd "/Users/fukui/Dropbox (Personal)/teaching/704_2025/704_Julia/Topic7/Stata"
use "JSTdatasetR6.dta",clear
bysort iso (year): gen inflation = 100*(log(cpi[_n+1]) - log(cpi))
bysort iso (year): gen GDP_growth = 100*(log(rgdpmad[_n]) - log(rgdpmad[_n-1]))
gen tloan_gdp = tloan/gdp

gen lgdp = log(rgdpmad)
gen lgdp_barro = log(rgdpbarro)

local country "USA"
if "`country'" == "JPN"{
	local yl "ylab(0(5)10)"
}
else if "`country'" == "USA"{
	local yl "ylab(0(5)20)"
}

tw (line stir year,lw(0.8) ) ///
(line inflation year, lw(0.8) lp(dash)) ///
if iso == "`country'" & inrange(year,1980,2019) ///
, legend(order(1 "Nominal Interest Rate" 2 "Inflation") rows(2) position(2) ring(0)) ytitle("p.p.") xtitle("Year")  name(ipi,replace) ///
title("Interest Rates and Inflation") `yl'

tw (line tloan_gdp year, lw(0.8)) ///
if iso == "`country'" & inrange(year,1980,2019) ///
, title("Debt to GDP Ratio") xtitle("Year") ytitle("") name(dgdp,replace)

tw (line lgdp year, lw(0.8)) ///
if iso == "`country'" & inrange(year,1980,2019) ///
, title("log Real GDP per capita") xtitle("Year") ytitle("") name(lgdp,replace)

graph combine dgdp ipi lgdp, rows(1) xsize(12) ysize(5)
graph export "./figure/`country'.pdf",replace

gen tloan_gdp = tloan/gdp
gen thh_gdp = thh/gdp
gen tbus_gdp = tbus/gdp


tw (line tbus_gdp year)  (line thh_gdp year) if iso== "JPN" & inrange(year,1980,2019)


merge m:1 countrycode year using ./original_data/pwt1001.dta, nogen keep(1 3)

tw line tloan_gdp year if countrycode== "JPN"

tw (line thh_gdp year, lw(0.8)) (line tbus_gdp year, lw(0.8))  if countrycode== "USA" & !missing(thh_gdp) ///
, legend(order(1 "Household debt" 2 "Firm debt") rows(1) position(6)) ytitle("% of GDP") xtitle("Year") 
graph export ./figure/hh_debt.pdf

tw ( scatter cwtfp thh_gdp if year == 2014 & thh_gdp!=., mlabel(countrycode) ) ///
	( lfit cwtfp thh_gdp if year == 2014 & thh_gdp!=.) ///
, legend(position(6)) name(tfp,replace)
