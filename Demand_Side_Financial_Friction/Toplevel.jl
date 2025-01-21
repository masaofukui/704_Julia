using QuantEcon
using Parameters
using Plots
using Interpolations
using SparseArrays
using LinearAlgebra
using BenchmarkTools
include("functions_sub.jl")
include("functions_Bellman_iteration.jl")

fig_save = 1

function set_parameters(;beta = nothing, amin = -2)
    # income process
    yg = [0.4; 1.0];
    Ny = length(yg);
    job_finding = 1-exp(-0.4*3);
    job_losing = 1-exp(-0.02*3);
    ytran = [1-job_losing job_losing; job_finding 1-job_finding];

    # compute stationary distribution of income process
    #mc = MarkovChain(ytran)
    #ss = stationary_distributions(mc)[1]

    Na = 100;

    # grid for assets
    amin = -2;
    amax = 20;
    ag = range(amin,amax,length=Na);

    # risk aversion for utility function
    gamma = 1.0;

    return (
        yg = yg,
        ag = ag,
        Ny = Ny,
        Na = Na,
        beta = beta,
        amin = amin,
        amax = amax,
        tol = 1e-4,
        gamma = gamma,
        ytran = ytran,
    )
end

param = set_parameters()


beta = 0.98
r = 0.02;
@assert beta*(1+r) < 1
Bellman_result = solve_policy_EGM(param, beta, r)
@unpack c_pol, a_pol = Bellman_result

@unpack ag = param
plot(ag,a_pol)

ss_distribution = solve_ss_distribution(param,Bellman_result)

plot(ag,ss_distribution)

@btime solve_ss_distribution(param,Bellman_result)
