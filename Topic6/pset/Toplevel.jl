using QuantEcon
using Parameters
using Plots
using Plots.Measures
include("functions_model_objects.jl")
include("functions_solve_eqm.jl")

fig_save = true
default_colors = palette(:auto);

function set_parameters(;rho_z = 0.85)

    beta = 0.96;
    sigma_z = 0.56;
    sig = sqrt(1-rho_z.^2).*sigma_z;
    mu = 0.0;
    alph = 1/3;

    Nz = 1000;

    Markov_Chain_rouwenhorst = tauchen(Nz,rho_z,sig, mu,6.0)
    lzg = Markov_Chain_rouwenhorst.state_values
    zg = exp.(lzg)
    ztran = Markov_Chain_rouwenhorst.p
    L = 1.0;
    delta = 0.05;
    lambda = 1.5;
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
        tol = 1e-8,
        alph = alph,
        L = L,
        delta = delta,
        lambda = lambda
    )
end

param = set_parameters(rho_z = 0.85)
rho_vec = range(0.5,0.99,length=20)
ss_result = ss_comparative_statics(param, rho_vec)


default(; titlefontfamily = "Computer Modern",
    xguidefontfamily = "Computer Modern",
    yguidefontfamily = "Computer Modern", 
    legendfontfamily = "Computer Modern",
    titlefontsize = 20,
    xguidefontsize = 12,
    legendfontsize = 12,
    yguidefontsize = 12,
    xgrid = :none)


plt_Z = plot(rho_vec,ss_result.Z_vec,linewidth = 5,label=:none)
xlabel!("Autocorrelation, \$\\rho_z\$")
title!("Total Factor Productivity, \$Z\$")

plt_Y = plot(rho_vec,ss_result.Y_vec, label = "Output",linewidth = 5,legend=:none)
xlabel!("Autocorrelation, \$\\rho_z\$")
title!("Output, \$Y\$")


plot(plt_Z,plt_Y,layout=(1,2),size=(1000,400))
plot!(margin=2mm)
if fig_save
    savefig("./Topic6/pset/figure/Z_Y_effect_of_rho_z.pdf")
end



T = 100;
omega_init = zeros(param.Nz);
omega_init[1] = 1.0;

plt_Y = plot();
plt_Z = plot();

rho_vec = [0.5,0.9,0.96];
for rho in rho_vec
    param = set_parameters(rho_z = rho)
    eqm_result = solve_eqm_path(param,omega_init, T; steady_state = false)
    plot!(plt_Y ,1:T,eqm_result.Y_path,label="\$\\rho_z = $rho\$",linewidth=5)
    xlabel!("Time")
    title!("Output, \$Y\$")

    plot!(plt_Z, 1:T,eqm_result.Z_path,label="\$\\rho_z = $rho\$",linewidth=5)
    xlabel!("Time")
    title!("Total Factor Productivity, \$Z\$")
end
plot(plt_Y,plt_Z,layout=(1,2),size=(1000,400))
plot!(margin=5mm)
if fig_save
    savefig("./Topic6/pset/figure/Y_Z_effect_of_rho_z.pdf")
end

