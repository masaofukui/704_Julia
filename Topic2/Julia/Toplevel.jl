using Debugger
using SparseArrays
using Plots
using DelimitedFiles
using DataFrames
using LinearAlgebra
using JLD
using Optim
using NLsolve
using Statistics
using Printf
using LaTeXStrings
using Plots.PlotMeasures
using Roots

colplot = palette(:Oranges_4);
colplot_green = palette(:Greens_4);
colplot_blue = palette(:Blues_3);


function set_parameters(;
        gamma = 0.75,
        alph = 0.75,
        s = 0.02,
        beta = (0.96).^(1/12),
        b = 0.4
    )

    P = Dict{String,Any}()
    m = 1
    P["m"] = m;
    P["alph"] = alph;
    
    uss = 0.05
    zss =1;
    lambdaU = s.*(1-uss)./uss
    thetass = (lambdaU/m).^(1/(1-alph))
    wssfun(c) = (1-gamma)*b + gamma*(zss + thetass*c)

    solvefun(c) = c - beta./(1-beta.*(1-s)).*qfun(P,thetass).*(zss - wssfun(c));
    c = fzero(solvefun,0,10)
    P["c"] = c;
    P["xi"] = 0.5;
    P["gamma"] = gamma;
    P["s"] = s;
    P["b"] = b;
    P["beta"] = beta;
    P["rho"] = 0.97327566
    P["uss"] = uss
    P["thetass"] = thetass;
    vss = uss*thetass;
    P["vss"] = vss
    P["omega"] = qfun(P,thetass)*vss;
    P["z shock"] = 1;
    return P
end


P = set_parameters()

include("./functions_sub.jl")
include("./Sequence_Space_Solver.jl")
include("./DMP_Model.jl")
include("./compute_steady_state.jl")

T = 500;

Shock_Var_in = "z shock"
Shock_Path = exp.(zeros(T))
Shock_Path[1] =  exp.(0.01)
for t = 2:T
    Shock_Path[t] = exp(P["rho"]*log(Shock_Path[t-1]))
end
ss0 = 0.1*ones(length(references))
ss0 = steady_state(ss0,references,P,varargin_eq,"z shock",ones(T))
ss1 = copy(ss0)
xfull = system_solve(ss0,ss1,references,P,varargin_eq,T,
Shock_Var_in,Shock_Path);

d = construct_IRF(xfull,n)

plot(1:T,d["ud"])




data, header = readdlm("./data/labor_productivity_bk.csv",',', header= true);
df =  DataFrame(data, vec(header));

lp_shock = data[:,3]
Tdata = length(lp_shock)
upath = zeros(Tdata)
for t = eachindex(lp_shock)
    endt = t:min(t+T-1,Tdata)
    IRF = d["ud"][1:length(endt)] .- P["uss"]
    upath[endt] +=  IRF./0.01.*lp_shock[t]
end

upath = upath .+ P["uss"]
plot(1:Tdata,upath)
open("./data/u_sim.csv", "w") do io
    writedlm(io, upath)
end



Shock_Var_in = "s"
Shock_Path = P["s"]*exp.(zeros(T))
Shock_Path[1] =  P["s"] + 0.01

ss0 = steady_state(ss0,references,P,varargin_eq,"z shock",ones(T))
ss1 = copy(ss0)
xfull = system_solve(ss0,ss1,references,P,varargin_eq,T,
Shock_Var_in,Shock_Path);
d = construct_IRF(xfull,n)



Tlim = 30;
lwset = 4
plt_u = plot(1:T,d["ud"], lw=lwset)
plot!(1:T,P["uss"]*ones(T), lw=lwset, color=:grey, linestyle=:dash )
xlims!((0, Tlim))
plot!(legend=:false)
plot!(xlabel="Months")
plot!(title="Unemployment rate, "*L"u")


plt_s = plot(1:T,Shock_Path, lw=lwset)
xlims!((0, Tlim))
plot!(1:T,P["s"]*ones(T), lw=lwset, color=:grey, linestyle=:dash )
plot!(legend=:false)
plot!(xlabel="Months")
plot!(title="Separation rate, "*L"s")


plt_v = plot(1:T,d["vd"], lw=lwset)
plot!(1:T,P["vss"]*ones(T), lw=lwset, color=:grey, linestyle=:dash )
xlims!((0, Tlim))
plot!(legend=:false)
plot!(xlabel="Months")
plot!(title="Vacancy, "*L" v")


plt_theta = plot(1:T,d["theta"], lw=lwset)
plot!(1:T,P["thetass"]*ones(T), lw=lwset, color=:grey, linestyle=:dash )
xlims!((0, Tlim))
plot!(legend=:false)
plot!(xlabel="Months")
plot!(title="Market tightness, "*L"\theta")
ylims!((0,0.05))

plot(plt_s,plt_u,plt_v,plt_theta,size=(800,500))
plot!(bottom_margin=4mm)
fig_name_save = "./figure/IRF_s.pdf"
savefig(fig_name_save)

##########################################################################################
dt = 0.01;
Tex = Int(1/dt)
up_ck = P["uss"]*ones(T*Tex)
up_ck[1] = P["uss"] + Shock_Path[1]*(1-P["uss"])
vp_ck = P["vss"]*ones(T*Tex)
for t = 2:(T*Tex)
    theta_temp = vp_ck[t-1]/up_ck[t-1]
    up_ck[t] = up_ck[t-1] + dt*(P["s"]*(1-up_ck[t-1]) - ffun(P,theta_temp)*up_ck[t-1]);
    vp_ck[t] = vp_ck[t-1] + dt*(P["omega"] - qfun(P,theta_temp)*vp_ck[t-1]);
end
tgrid = dt*(1:(T*Tex))
plot(tgrid,vp_ck)
plot(tgrid,up_ck)

theta_ck = vp_ck./up_ck


Tlim = T;
lwset = 4
plt_u = plot(tgrid,up_ck, lw=lwset)
plot!(1:T,P["uss"]*ones(T), lw=lwset, color=:grey, linestyle=:dash )
xlims!((0, Tlim))
plot!(legend=:false)
plot!(xlabel="Months")
plot!(title="Unemployment rate, "*L"u")


plt_s = plot(1:T,Shock_Path, lw=lwset)
xlims!((0, Tlim))
plot!(1:T,P["s"]*ones(T), lw=lwset, color=:grey, linestyle=:dash )
plot!(legend=:false)
plot!(xlabel="Months")
plot!(title="Separation rate, "*L"s")


plt_v = plot(tgrid,vp_ck, lw=lwset)
plot!(1:T,P["vss"]*ones(T), lw=lwset, color=:grey, linestyle=:dash )
xlims!((0, Tlim))
plot!(legend=:false)
plot!(xlabel="Months")
plot!(title="Vacancy, "*L" v")


plt_theta = plot(tgrid,theta_ck, lw=lwset)
plot!(1:T,P["thetass"]*ones(T), lw=lwset, color=:grey, linestyle=:dash )
xlims!((0, Tlim))
plot!(legend=:false)
plot!(xlabel="Months")
plot!(title="Market tightness, "*L"\theta")
ylims!((0,0.05))

plot(plt_s,plt_u,plt_v,plt_theta,size=(800,500))
plot!(bottom_margin=4mm)
fig_name_save = "./figure/IRF_s_ck.pdf"
savefig(fig_name_save)
