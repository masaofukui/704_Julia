using QuantEcon
using Parameters
using Plots
using Interpolations
using SparseArrays
using LinearAlgebra
using Roots
include("functions_sub.jl")
include("functions_Bellman_iteration.jl")
include("functions_solve_distribution.jl")
function set_parameters(;r = nothing, beta = nothing, phi = 2.0, Y = 1.0,a_range=5.0)
    # income process
    yg = [0.4; 1.0];
    Ny = length(yg);
    job_finding = 1-exp(-0.4*3);
    job_losing = 1-exp(-0.02*3);
    ytran = [1-job_finding job_finding ; job_losing 1-job_losing ];
    
    # make sure ytran is a transition matrix
    @assert sum(ytran,dims=2) ≈ ones(Ny)
    # compute stationary distribution of income process
    #mc = MarkovChain(ytran)
    #ss = stationary_distributions(mc)[1]

    Na = 100;
    T = 300; 

    # grid for assets
    amin = -phi;
    amax = amin + a_range;
    ag = range(amin,amax,length=Na);
    ag_plot = ag .- phi;

    # risk aversion for utility function
    gamma = 1.0;
    # step size for numerical differentiation
    dx = 0.01; 
    return (
        yg = yg,
        ag = ag,
        ag_plot = ag_plot,
        Ny = Ny,
        Na = Na,
        beta = beta,
        amin = amin,
        amax = amax,
        tol = 1e-4,
        gamma = gamma,
        ytran = ytran,
        T,
        dx, 
        phi = phi,
        r = r,
        Y = Y,
    )
end


########################################################
# Bellman Iteration
########################################################
beta = 0.96
r = 0.02
param = set_parameters(r = r, beta = beta, phi = 2.0,  Y = 1.0,a_range=5.0)
Bellman_result = solve_policy_EGM(param, beta, r)

plot(param.ag,Bellman_result.c_pol,
    lw = 5,
    xlabel = "Assets",
    ylabel = "Consumption",
    title = "Consumption Function",
    label = ["low income" "high income"]
)
plot(param.ag,Bellman_result.a_pol,
    lw = 5,
    xlabel = "Assets",
    ylabel = "Asset next period",
    title = "Asset Policy Function",
    label = ["low income" "high income"]
)

########################################################
# Stationary Distribution
########################################################
ss_distribution = solve_ss_distribution(param,Bellman_result)

plot(param.ag,ss_distribution,
    lw = 5,
    xlabel = "Assets",
    ylabel = "Stationary Distribution",
    title = "Stationary Distribution"
)

########################################################
# Asset Demand
########################################################

function compute_asset_demand(param,beta,r)
    @unpack ag = param
    Bellman_result = solve_policy_EGM(param, beta, r)
    ss_distribution = solve_ss_distribution(param,Bellman_result)
    A = sum(ss_distribution.*ag)
    return A
end

solve_r = x -> compute_asset_demand(param,beta,x) - 0.0
r_ss_eqm = find_zero(solve_r, (-0.1, 1/beta-1-1e-5))



