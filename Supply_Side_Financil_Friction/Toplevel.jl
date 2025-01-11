using QuantEcon
using Parameters
using Plots
include("functions_model_objects.jl")
include("functions_solve_eqm.jl")



function set_parameters()

    beta = 0.96;
    rho_z = 0.99;
    sigma_z = 2;
    sig = (1-rho_z)*sigma_z;
    mu = 0;

    Nz = 100;

    Markov_Chain_rouwenhorst = rouwenhorst(Nz,rho_z,sig, mu)
    lzg = Markov_Chain_rouwenhorst.state_values
    zg = exp.(lzg)
    ztran = Markov_Chain_rouwenhorst.p

    #=
    ss = stationary_distributions(Markov_Chain_rouwenhorst)[1]
    plot(zg,ss)
    =#
    return (
        beta = beta,
        Nz = Nz,
        zg = zg,
        ztran = ztran,
        tol = 1e-6,
    )
end

param = set_parameters()
T = 1000;
lambda_path = 3*ones(T);
omega_init = ones(param.Nz)/param.Nz;
eqm_result = solve_eqm_path(param,omega_init,lambda_path, T)
eqm_result.Tend

eqm_result.r_index_path

plot(param.zg,eqm_result.omega_path[:,end], label = "r")
plot!(param.zg,eqm_result.omega_path[:,end-1], label = "r")
maximum(abs.(eqm_result.omega_path[:,100] - eqm_result.omega_path[:,100-1]))
plot(1:T,eqm_result.Z_path, label = "Z")