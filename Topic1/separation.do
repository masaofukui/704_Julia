* append your graph command to this file: e.g.
* tsline timeseriesvar, xlabel(,format(%tm))  legend(off) graphregion(color(white)) legend(order(6 1 "Recession"))
local yms = ym(1980,1)
local ymd = ym(2020,1)
twoway function y=3,range(240 246) recast(area) color(gs12) base(1) || /// 
function y=3,range(258 274) recast(area) color(gs12) base(1) || /// 
function y=3,range(366 374) recast(area) color(gs12) base(1) || /// 
function y=3,range(494 502) recast(area) color(gs12) base(1) || /// 
function y=3,range(575 593) recast(area) color(gs12) base(1) || /// 
(tsline separation if tin(1980m1,2019m12) , lw(0.5)  color(maroon)), ///
xlabel(`yms'(120)`ymd',format(%tmCCYY))   legend(off)  ///
graphregion(color(white)) legend(order(6 1 "Recession")) ///
ytitle("%") title("Separation rate, s") xsize(8)

graph export ./figure/separation.pdf,replace

local yms = ym(1980,1)
local ymd = ym(2020,1)
twoway function y=45,range(240 246) recast(area) color(gs12) base(15) || /// 
function y=45,range(258 274) recast(area) color(gs12) base(15) || /// 
function y=45,range(366 374) recast(area) color(gs12) base(15) || /// 
function y=45,range(494 502) recast(area) color(gs12) base(15) || /// 
function y=45,range(575 593) recast(area) color(gs12) base(15) || /// 
(tsline job_finding if tin(1980m1,2019m12) ,  lw(0.5) color(maroon)), ///
xlabel(`yms'(120)`ymd',format(%tmCCYY))   legend(off)  ///
graphregion(color(white)) legend(order(6 1 "Recession")) ///
ytitle("%") title("Job-finding rate, f") xsize(8) ///
ylabel(20(10)40)
graph export ./figure/job_finding.pdf,replace
