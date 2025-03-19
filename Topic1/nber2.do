* append your graph command to this file: e.g.
* tsline timeseriesvar, xlabel(,format(%tm))  legend(off) graphregion(color(white)) color(maroon) legend(order(6 1 "Recession"))
local yms = ym(1980,1)
local ymd = ym(2020,1)
twoway function y=.6181199884414673,range(240 246) recast(area) color(gs14) base(-.8) || /// 
function y=.6181199884414673,range(258 274) recast(area) color(gs14) base(-.8) || /// 
function y=.6181199884414673,range(366 374) recast(area) color(gs14) base(-.8) || /// 
function y=.6181199884414673,range(494 502) recast(area) color(gs14) base(-.8) || /// 
function y=.6181199884414673,range(575 593) recast(area) color(gs14) base(-.8) || /// 
( tsline de if tin(1980m1,2019m12) , color(maroon) lw(0.5)), xlabel(`yms'(120)`ymd',format(%tmCCYY))  legend(off) graphregion(color(white))  ///
 title("Net Inflow between U and E")  ylabel(-0.5(0.5)0.5) 
graph export ./figure/net_u_e.pdf,replace


local yms = ym(1980,1)
local ymd = ym(2020,1)
twoway function y=.6181199884414673,range(240 246) recast(area) color(gs14) base(-.8) || /// 
function y=.6181199884414673,range(258 274) recast(area) color(gs14) base(-.8) || /// 
function y=.6181199884414673,range(366 374) recast(area) color(gs14) base(-.8) || /// 
function y=.6181199884414673,range(494 502) recast(area) color(gs14) base(-.8) || /// 
function y=.6181199884414673,range(575 593) recast(area) color(gs14) base(-.8) || /// 
(tsline dn if tin(1980m1,2019m12) , lw(0.5) color(maroon)), xlabel(`yms'(120)`ymd',format(%tmCCYY))  legend(off) graphregion(color(white))  ///
 title("Net Inflow between U and N") ylabel(-0.5(0.5)0.5) 
 graph export ./figure/net_n_e.pdf,replace
