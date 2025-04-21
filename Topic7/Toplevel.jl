using QuantEcon
using Parameters
using Plots
using Interpolations
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using LaTeXStrings
using ForwardDiff
include("functions_sub.jl")
include("functions_Bellman_iteration_transorm_phi.jl")
include("functions_distribution_iteration.jl") 
include("functions_ssj.jl")
include("functions_plot.jl")
include("functions_simulation.jl")
include("functions_GE_steady_state.jl")
fig_save = 1
default_colors = palette(:auto);

function set_parameters(;r = nothing, beta = nothing, phi = 2.0, transform_phi = false, Y = 1.0)
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
    if transform_phi
        amin = 0;
    else
        amin = -phi;
    end
    amax = amin + 5.0;
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

fig_save = 1
beta = 0.978
r = 0.02;
param = set_parameters(r=r,beta=beta,transform_phi = true, phi=2.0)
param_changed = 
    (beta = param.beta, 
    r = r, 
    phi = param.phi,
    Y = param.Y)

########################################################
# Baseline results
########################################################
fig_name = "beta_"*string(beta)*"_r_"*string(r);
wrapper_PE_policy_plot(param,param_changed,fig_name = fig_name,fig_save = fig_save)
Bellman_result = solve_policy_EGM(param, param_changed)
plt_all = run_simulation(param, Bellman_result; seed_num = 1231,fig_save = fig_save)
ss_distribution = solve_ss_distribution(param,Bellman_result)

plot_distribution(param, Bellman_result, ss_distribution,fig_save = fig_save, fig_name = fig_name)


########################################################
# Comparative Statics w.r.t. phi
########################################################
phi_vec = [2.0, 1.5]
A_vec, Bellman_list, ss_distribution_list = compute_asset_demand_vec(param,phi_vec,var_change = "phi")
plot_distribution_compare(param, Bellman_list, ss_distribution_list,phi_vec,fig_save = fig_save,fig_name = "phi")


########################################################
# Comparative Statics w.r.t. r
########################################################
r_vec = [0.02; 0.021]
A_vec, Bellman_list, ss_distribution_list = compute_asset_demand_vec(param,r_vec,var_change = "r")
plot_distribution_compare(param, Bellman_list, ss_distribution_list,phi_vec,fig_save = fig_save,fig_name = "r")

########################################################
# Comparative Statics w.r.t. Y
########################################################
Y_vec = [1.0; 1.1]
A_vec, Bellman_list, ss_distribution_list = compute_asset_demand_vec(param,Y_vec,var_change = "Y")
plot_distribution_compare(param, Bellman_list, ss_distribution_list,Y_vec,fig_save = fig_save,fig_name = "Y",var_change = "Y")

########################################################
# Plot the Asset Demand Function
########################################################
beta = 0.97
r = 0.02;
phi = 0.8
param = set_parameters(r=r,beta=beta,transform_phi = true, phi=phi)

r_vec = range(0.0,1/param.beta-1 - 1e-4,length=40)
A_vec, Bellman_list, ss_distribution_list = compute_asset_demand_vec(param,r_vec,var_change = "r")
A_supply = zeros(length(r_vec))
plot(r_vec,A_vec,
    lw=5,
    label="Asset Demand, \$A(r)\$",
    xlabel="Interest rate, \$r\$"
)
plot!(r_vec,A_supply,
    lw=5,
    label="Asset Supply",
    linestyle=:dash
)



########################################################
# Calibrate beta
########################################################
r_target = 0.02
phi = 0.8
beta_calibrated = calibrate_beta(param,r_target)
param_calibrated = set_parameters(r=r_target,beta=beta_calibrated,transform_phi = true, phi=phi)

########################################################
# Steady State Counterfactual
########################################################
@time compute_r(param_calibrated)

#trying out the functions 
c_ss, a_ss, ss_distribution = get_ss_objects(param, beta, r)
Lambda =  construct_transition_matrix(param, a_ss)
E = get_expectations(param, Lambda, beta, r, a_ss)
curlyY, curlyD = get_curly_Y_D(param, beta, r, c_ss, a_ss, ss_distribution)
 J =  get_SSJ(param, curlyY, curlyD, E)[1]
 F =  get_SSJ(param, curlyY, curlyD, E)[2]
