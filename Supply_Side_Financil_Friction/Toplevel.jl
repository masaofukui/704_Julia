using QuantEcon
using Parameters
using Plots
include("functions_model_objects.jl")



function set_parameters()

    beta = 0.96;
    rho_z = 0.99;
    sigma_z = 2;
    sig = (1-rho_z)*sigma_z;
    mu = 0;

    Nz = 500;

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
        lambda = lambda,
        Delta_z,
        tol = 1e-6,
    )
end

function solve_eqm_path(param,omega_init,lambda_path, T)
    @assert length(lambda_path) == T
    @unpack Nz,tol = param

    omega_path = zeros(Nz,T)
    omega_path[:,1] = omega_init;
    r_path = zeros(T)
    r_index_path = zeros(T)
    Z_path = zeros(T)
    Kgrowth_path = zeros(T)

    Tend = []
    for t = 1:T
        omega = omega_path[:,t];
        lambda = lambda_path[t];
        r, r_index = solve_for_r(param,omega,lambda)
        piz = compute_piz(param,r,lambda)
        Z = compute_Z(param,omega,r_index);
        Kgrowth = compute_Kgrowth(param,omega,piz)
        omega_next = compute_omega_next(param,omega,piz,Kgrowth)

        r_path[t] = r;
        r_index_path[t] = r_index;
        Z_path[t] = Z;
        Kgrowth_path[t] = Kgrowth;
        if t < T
            omega_path[:,t+1] = omega_next;
            if maximum(abs.(omega_path[:,t+1] .- omega_path[t])) < tol
                println("Steady state reached at t = $t")
                Tend = copy(t);
                break
            else
                Tend = copy(T);
            end
        end        
    end

    return (
        omega_path = omega_path[:,1:Tend],
        r_path = r_path[1:Tend],
        r_index_path = r_index_path[1:Tend],
        Z_path = Z_path[1:Tend],
        Kgrowth_path = Kgrowth_path[1:Tend],
        Tend = Tend
    )
end

param = set_parameters()
T = 1000;
lambda_path = 3*ones(T);
omega_init = ones(param.Nz)/param.Nz;
eqm_result = solve_eqm_path(param,omega_init,lambda_path, T)
eqm_result.Tend

plot(param.zg,eqm_result.omega_path, label = "r")

plot(param.zg,eqm_result.omega_path[:,end], label = "r")

plot(1:T,eqm_result.Z_path, label = "Z")