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


function set_parameters(;r = nothing, beta = nothing, phi = 2.0, Y = 1.0,a_range=1e7,
    return_heterogeneity = false,
    nonhomothetic = false,
    scale_dependent = false,
    )
    # income process

    # Define continuous Pareto distribution
    zeta_y = 2.2  # shape
    ymin = 0.1 # scale (minimum value)
    pareto_y = Pareto(zeta_y,ymin)

    # Define grid over probabilities
    Ny = 50  # number of points
    ymax = quantile.(Ref(pareto_y), 0.99999)


    # Get quantiles
    lyg = range(log(ymin),log(ymax),length=Ny)
    yg = exp.(lyg)
    yg = yg.*Y;
    cdf_y = cdf(pareto_y,yg)
    pdf_y = diff(cdf_y)
    pdf_y_res = 1-sum(pdf_y)
    pdf_y = [pdf_y; pdf_y_res]

    if return_heterogeneity
        zg = [0.0,2.0];
        z_probs = [0.5,0.5];
        death_prob = 0.00;
    else
        zg = 1.0;
        z_probs = 1.0;
        death_prob = 0.00;
    end
    Nz = length(zg)

   
    p_draw = 0.2;
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
    gamma = 4.0;
    # wealth utility function
    if nonhomothetic
        nu = 0.6;
        kappa = 1e-9;
        death_prob = 0.05;
    else
        nu = 0.0;
        kappa = 0.0;
    end
    if scale_dependent
        eta = 0.15;
        psi = 1e-12;
    else
        eta = 0.0;
        psi = 1.0;
    end
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
        nu = nu,
        kappa = kappa,
        eta = eta,
        psi = psi,
    )
end

beta = 0.978
r = 0.02;

########################################################
# Baseline results
########################################################
param = set_parameters(r=r,beta=beta, phi=0.0)
fig_name = "beta_"*string(beta)*"_r_"*string(r);
Bellman_result = solve_policy_EGM_z(param, beta,r)
ss_distribution = solve_ss_distribution(param,Bellman_result)


plot_power_law_wrapper(param,Bellman_result,ss_distribution,fig_save = fig_save, fig_name = "canonical_all_distributions")

########################################################
# Return heterogeneity
########################################################

param = set_parameters(r=r,beta=beta, phi=0.0,return_heterogeneity = true)
Bellman_result = solve_policy_EGM_z(param, beta,r)
plot(param.ag,Bellman_result.c_pol[:,:,1],label=:none)
ss_distribution = solve_ss_distribution(param,Bellman_result)
plot(log.(param.ag),(vec(sum(ss_distribution,dims=[2,3]))),label="distribution")

plot_power_law_wrapper(param,Bellman_result,ss_distribution,fig_save = fig_save, fig_name = "return_heterogeneity")


########################################################
# Nonhomothetic and scale dependent
########################################################
beta = 0.3
r = 0.02;
param = set_parameters(r=r,beta=beta, phi=0.0,return_heterogeneity = true,nonhomothetic = true,scale_dependent = true)
Bellman_result = solve_policy_EGM_z(param, beta,r)
plot(param.ag,Bellman_result.c_pol[:,:,1],label=:none)

ss_distribution = solve_ss_distribution(param,Bellman_result)
plot(log.(param.ag),(vec(sum(ss_distribution,dims=[2,3]))),label="distribution")
plot_power_law_wrapper(param,Bellman_result,ss_distribution,fig_save = fig_save, fig_name = "return_heterogeneity_nonhomothetic_scale_dependent")
