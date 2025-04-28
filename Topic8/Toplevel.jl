using QuantEcon
using Parameters
using Plots
using Interpolations
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using LaTeXStrings
using ForwardDiff
using Distributions
using Plots.Measures
include("functions_sub.jl")
include("functions_Bellman_iteration_z.jl")
include("functions_distribution_iteration.jl") 
include("functions_plot.jl")
fig_save = 1
default_colors = palette(:auto);
fig_dir = "./figure/Topic8/"


function set_parameters(;r = nothing, beta = nothing, phi = 2.0, Y = 1.0,a_range=1e5,
    return_heterogeneity = false)
    # income process

    # Define continuous Pareto distribution
    zeta_y = 2.2  # shape
    ymin = 0.1 # scale (minimum value)
    pareto_y = Pareto(zeta_y,ymin)

    # Define grid over probabilities
    Ny = 100  # number of points
    probs = range(0, stop=0.99999, length=Ny)  # avoid 0 and 1 exactly

    # Get quantiles
    yg = quantile.(Ref(pareto_y), probs)
    yg = yg.*Y;
    cdf_y = cdf(pareto_y,yg)
    pdf_y = diff(cdf_y)
    pdf_y_res = 1-sum(pdf_y)
    pdf_y = [pdf_y; pdf_y_res]

    if return_heterogeneity
        zg = [0.0,5];
        z_probs = [0.5,0.5];
        death_prob = 0.05;
    else
        zg = 1.0;
        z_probs = 1.0;
        death_prob = 0.0;
    end
    Nz = length(zg)

   
    p_draw = 0.1;
    ytran = p_draw*ones(Ny,1)*pdf_y' + (1-p_draw)*I;
    
    # make sure ytran is a transition matrix
    @assert sum(ytran,dims=2) ≈ ones(Ny)
    # compute stationary distribution of income process
    #mc = MarkovChain(ytran)
    #ss = stationary_distributions(mc)[1]

    Na = 200;

    # grid for assets

    amin = -phi;
    amax = amin + a_range;
    ag = discretize_assets(amin,amax,Na);
    ag_plot = ag .- phi;

    # risk aversion for utility function
    gamma = 1.0;
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
        phi = phi,
        r = r,
        Y = Y,
        zg = zg,
        z_probs = z_probs,
        death_prob = death_prob,
        Nz = Nz,
    )
end

beta = 0.978
r = 0.02;
param = set_parameters(r=r,beta=beta, phi=0.0)

########################################################
# Baseline results
########################################################
fig_name = "beta_"*string(beta)*"_r_"*string(r);
Bellman_result = solve_policy_EGM_z(param, beta,r)
ss_distribution = solve_ss_distribution(param,Bellman_result)


plot_power_law_wrapper(param,Bellman_result,ss_distribution,fig_save = fig_save, fig_name = "canonical_all_distributions")

########################################################
# Return heterogeneity
########################################################
param = set_parameters(r=r,beta=beta, phi=0.0,return_heterogeneity = true)
Bellman_result = solve_policy_EGM_z(param, beta,r)
ss_distribution = solve_ss_distribution(param,Bellman_result)

plot_power_law_wrapper(param,Bellman_result,ss_distribution,fig_save = fig_save, fig_name = "canonical_all_distributions")
