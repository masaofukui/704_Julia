using QuantEcon
using Parameters
using Plots
include("functions_model_objects.jl")
include("functions_solve_eqm.jl")

fig_save = true
default_colors = palette(:auto);

function set_parameters()

    beta = 0.96;
    rho_z = 0.85;
    sigma_z = 0.56;
    sig = sqrt(1-rho_z.^2).*sigma_z;
    mu = 0;

    Nz = 500;

    Markov_Chain_rouwenhorst = tauchen(Nz,rho_z,sig, mu,4.0)
    lzg = Markov_Chain_rouwenhorst.state_values
    zg = exp.(lzg)
    ztran = Markov_Chain_rouwenhorst.p

    #=
    ss = stationary_distributions(Markov_Chain_rouwenhorst)[1]
    plot(zg,ss)
    Elogz = sum(log.(zg) .*ss)
    Var_logz = sum((log.(zg) .-Elogz).^2 .*ss)
    std_logz = sqrt(Var_logz)
    =#
    return (
        beta = beta,
        Nz = Nz,
        zg = zg,
        ztran = ztran,
        tol = 1e-4,
    )
end

param = set_parameters()
lambda_vec = 1:0.25:3;
ss_result = ss_comparative_statics(param, lambda_vec)
direct_effect_result = ss_comparative_statics(param, lambda_vec,T=1,omega_in = ss_result.omega_vec[:,1])

default(; titlefontfamily = "Computer Modern",
    xguidefontfamily = "Computer Modern",
    yguidefontfamily = "Computer Modern", 
    legendfontfamily = "Computer Modern",
    titlefontsize = 20,
    xguidefontsize = 12,
    legendfontsize = 12,
    yguidefontsize = 12,
    xgrid = :none)


plot(lambda_vec,ss_result.Z_vec, label = "Total",linewidth = 5,legend=:none)
xlabel!("Borrowing constraint parameter, \$\\lambda\$")
title!("Total Factor Productivity, Z")
plot!(xgrid=:none)

if fig_save
    savefig("./figure/Z_effect_of_lambda.pdf")
end
plot!(lambda_vec,direct_effect_result.Z_vec, label = "Direct Effect",linewidth = 5,linestyle=:dash,legend= :topleft)


plot(lambda_vec,ss_result.Z_vec, label = "Long-run Effect",linewidth = 0,legend=:topleft,color=default_colors[1],fillrange=minimum(direct_effect_result.Z_vec),fillalpha=1.0)
plot!(lambda_vec,direct_effect_result.Z_vec, label = "Impact Effect",linewidth = 0,color=default_colors[2],fillrange=minimum(direct_effect_result.Z_vec),fillalpha=1.0)
xlabel!("Borrowing constraint parameter, \$\\lambda\$")
title!("Total Factor Productivity, Z")
if fig_save
    savefig("./figure/Z_decomposition.pdf")
end


plot(param.zg,ss_result.omega_vec[:,1],linewidth=5,label="Wealth distribution, \$\\omega(z)\$")
plot!(param.zg,ss_result.kz_vec[:,1],linewidth=5,linestyle=:dash,label="Capital distribution, \$ k(z)\$")
xlabel!("Productivity, z")
title!("\$\\lambda\$ = $(lambda_vec[1])")
if fig_save
    savefig("./figure/autarky_wealth_capital_distribution.pdf")
end


ilambda = 3
plot(param.zg,ss_result.omega_vec[:,1],linewidth=5,label="Wealth distribution, \$\\omega(z)\$")
plot!(param.zg,direct_effect_result.kz_vec[:,ilambda],linewidth=5,linestyle=:dash,label="Capital distribution, \$ k(z)\$")
xlabel!("Productivity, z")
title!("\$\\lambda\$ = $(lambda_vec[ilambda])")
if fig_save
    savefig("./figure/autarky_wealth_capital_distribution_lambda_$(lambda_vec[ilambda]).pdf")
end


plot(param.zg,ss_result.omega_vec[:,1],linewidth=5,label=:none,alpha=0.2)
plot!(param.zg,direct_effect_result.kz_vec[:,ilambda],linewidth=5,linestyle=:dash,label=:none,alpha=0.0,color=default_colors[2])
plot!(param.zg,ss_result.omega_vec[:,ilambda],linewidth=5,label="Wealth distribution, \$\\omega(z)\$",color=default_colors[1])
xlabel!("Productivity, z")
title!("\$\\lambda\$ = $(lambda_vec[ilambda])")
if fig_save
    savefig("./figure/long_run_wealth_distribution_lambda_$(lambda_vec[ilambda]).pdf")
end
plot!(param.zg,direct_effect_result.kz_vec[:,ilambda],linewidth=5,linestyle=:dash,label=:none,alpha=0.2,color=default_colors[2])
plot!(param.zg,ss_result.kz_vec[:,ilambda],linewidth=5,linestyle=:dash,label="Capital distribution, \$ k(z)\$",color=default_colors[2])
if fig_save
    savefig("./figure/long_run_wealth_capital_distribution_lambda_$(lambda_vec[ilambda]).pdf")
end


##################################################################
# short-run effect
##################################################################
ilambda = 3
T = 30;
omega_init = ss_result.omega_vec[:,ilambda];
lambda_ss = lambda_vec[ilambda];
lambda_path = zeros(T);
lambda_path[1] = 1.3;
for t in 1:(T-1)
    lambda_path[t+1] = lambda_ss + 0.5*(lambda_path[t] - lambda_ss);
end
tpre = 3;
tg = -tpre:T;


eqm_result = solve_eqm_path(param,omega_init,lambda_path,T; steady_state = false);



lambda_plot = lambda_ss*ones(T+tpre+1);
lambda_plot[tpre+2:end] = lambda_path;
plot(tg,lambda_plot,linewidth=5,label=:none)
xlabel!("Year")
title!("Borrowing constraint, \$\\lambda\$")
plot!(size=(500,500))
if fig_save
    savefig("./figure/lambda_path.pdf")
end

Z_plot = [eqm_result.Z_path[end]*ones(tpre+1); eqm_result.Z_path];
plot(tg,Z_plot,label=:none,linewidth=5)
xlabel!("Year")
title!("Total Factor Productivity, \$Z\$")
plot!(size=(500,500))

if fig_save
    savefig("./figure/Z_effect_of_lambda_short_run.pdf")
end

ss_growth = eqm_result.Kgrowth_path[end];
lY_plot = zeros(T+tpre+1);
lK_plot = zeros(T+tpre+1);
lK_plot[1] = 0.0;
lY_no_shock = zeros(T+tpre+1);
lK_no_shock = zeros(T+tpre+1);
lK_no_shock[1] = 0.0;
for it in eachindex(tg)
    t = tg[it];
    if t < 0 
        lK_plot[it+1] = lK_plot[it] + log(ss_growth);
    elseif it < length(lK_plot)
        println("t = $t, Kgrowth = $(eqm_result.Kgrowth_path[t+1])")
        lK_plot[it+1] = lK_plot[it] + log(eqm_result.Kgrowth_path[t+1]);
    end
    if it < length(lK_no_shock)
        lK_no_shock[it+1] = lK_no_shock[it] + log(ss_growth);
    end
    lY_plot[it] = log(Z_plot[it]) + lK_plot[it]
    lY_no_shock[it] = log(Z_plot[1]) + lK_no_shock[it]
end

plot(tg,lY_plot .- lY_no_shock,label=:none,linewidth=5)
xlabel!("Year")
title!("Log Output, \$\\log(Y)\$")
plot!(size=(500,500))
if fig_save
    savefig("./figure/output_effect_of_lambda_short_run.pdf")
end










