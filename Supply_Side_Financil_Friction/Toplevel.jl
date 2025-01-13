using QuantEcon
using Parameters
using Plots
include("functions_model_objects.jl")
include("functions_solve_eqm.jl")

fig_save = 1

function set_parameters()

    beta = 0.96;
    rho_z = 0.95;
    sigma_z = 0.1;
    sig = sqrt(1-rho_z.^2).*sigma_z;
    mu = 0;

    Nz = 100;

    Markov_Chain_rouwenhorst = rouwenhorst(Nz,rho_z,sig, mu)
    lzg = Markov_Chain_rouwenhorst.state_values
    zg = exp.(lzg)
    ztran = Markov_Chain_rouwenhorst.p

    #=
    ss = stationary_distributions(Markov_Chain_rouwenhorst)[1]
    plot(zg,ss)
    Elogz = sum(log.(zg) .*ss)
    Varlogz = sum((log.(zg) .-Elogz).^2 .*ss)
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
lambda_vec = 1:0.5:10;
ss_result = ss_comparative_statics(param, lambda_vec)
direct_effect_result = ss_comparative_statics(param, lambda_vec,T=1,omega_in = ss_result.omega_vec[:,1])




plot(lambda_vec,ss_result.Z_vec, label = "Total",linewidth = 5,legend=:none)
xlabel!("Borrowing constraint parameter, \$\\lambda\$")
title!("Total Factor Productivity, Z")
plot!(xgrid=:none)
plot!(titlefontfamily = "Computer Modern",
    xguidefontfamily = "Computer Modern",
    yguidefontfamily = "Computer Modern",
    legendfontfamily = "Computer Modern",
    titlefontsize=20,xguidefontsize=12,legendfontsize=12,yguidefontsize=12)
if fig_save == 1
    savefig("./figure/Z_effect_of_lambda.pdf")
end
plot!(lambda_vec,direct_effect_result.Z_vec, label = "Direct Effect",linewidth = 5,linestyle=:dash,legend= :topleft)

if fig_save == 1
    savefig("./figure/Z_effect_of_lambda_direct.pdf")
end



plot(param.zg,ss_result.omega_vec[:,1],linewidth=5,label="\$\\lambda\$ = $(lambda_vec[1])")
plot!(param.zg,ss_result.omega_vec[:,end],linewidth=5,linestyle=:dash,label="\$\\lambda\$ = $(lambda_vec[end])")
xlabel!("Productivity, z")
title!("Wealth Share, \$ \\omega(z)\$")
plot!(xgrid=:none)
plot!(titlefontfamily = "Computer Modern",
    xguidefontfamily = "Computer Modern",
    yguidefontfamily = "Computer Modern",
    legendfontfamily = "Computer Modern",
    titlefontsize=20,xguidefontsize=12,legendfontsize=12,yguidefontsize=12)
if fig_save == 1
    savefig("./figure/omega_dist.pdf")
end
