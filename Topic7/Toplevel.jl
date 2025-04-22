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
include("functions_GE_transitions.jl")
include("functions_iMPC.jl")
fig_save = 1
default_colors = palette(:auto);
fig_dir = "./figure/Topic7/"

function set_parameters(;r = nothing, beta = nothing, phi = 2.0, transform_phi = false, Y = 1.0,a_range=5.0,higher_risk = false)
    # income process
    yg = [0.4; 1.0];
    Ny = length(yg);
    job_finding = 1-exp(-0.4*3);
    job_losing = 1-exp(-0.02*3);
    if higher_risk
        job_finding = 0.5;
        job_losing = 0.2;
    end
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
Bellman_result = solve_policy_EGM(param, param_changed = param_changed)
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
plot_asset_demand_function(param, r_vec, A_vec)

phi = 0.4
param_new = set_parameters(r=r,beta=beta,transform_phi = true, phi=phi)
r_vec_new = range(0.0,1/param_new.beta-1 - 1e-4,length=40)
A_vec_new, Bellman_list_new, ss_distribution_list_new = compute_asset_demand_vec(param_new,r_vec_new,var_change = "r")
plot_asset_demand_function(param, r_vec, A_vec_new; A_vec_old = A_vec,compare = true)

########################################################
# Calibrate beta
########################################################
r_target = 0.02
phi = 0.8
param = set_parameters(r=r_target,transform_phi = true, phi=phi)
beta_calibrated = calibrate_beta(param,r_target)
param_calibrated = set_parameters(r=r_target,beta=beta_calibrated,transform_phi = true, phi=phi)


########################################################
# Steady State Counterfactual: changes in phi
########################################################
phi_vec = range(0.8,0.2,length=10)
r_vec_counterfactual = zeros(length(phi_vec))
for (i,phi) in enumerate(phi_vec)
    param_in = set_parameters(r=r_target,beta=beta_calibrated,transform_phi = true, phi=phi)
    r_vec_counterfactual[i] = compute_r(param_in)
end

plt_r = plot(phi_vec,r_vec_counterfactual,
    lw=5,
    label=:none,
    xlabel="Borrowing Limit, \$\\phi\$",
    title="Interest Rate, \$r\$"
)
if fig_save == 1
    savefig(plt_r, fig_dir*"r_ss_counterfactual.pdf")
end

########################################################
# Steady State Counterfactual: changes in phi under demand determined equilibrium
########################################################
phi_vec = range(0.8,0.2,length=10)
Y_vec_counterfactual = zeros(length(phi_vec))
for (i,phi) in enumerate(phi_vec)
    param_in = set_parameters(r=r_target,beta=beta_calibrated,transform_phi = true, phi=phi)
    Y_vec_counterfactual[i] = compute_Y(param_in)
end

plt_Y = plot(phi_vec,Y_vec_counterfactual,
    lw=5,
    label=:none,
    xlabel="Borrowing Limit, \$\\phi\$",
    title="Output, \$Y\$",
    color = default_colors[1]
)
plot!(phi_vec,ones(length(phi_vec))*param.Y,
    lw=5,
    label=:none,
    alpha = 0.2,
    color = default_colors[1],
    ls = :dash
)
plt_r = plot(phi_vec,ones(length(phi_vec))*param.r,
    lw=5,
    color = default_colors[1],
    ls = :solid,
    label = "Fixed \$r\$"
)
plot!(phi_vec,r_vec_counterfactual,
    lw=5,
    xlabel="Borrowing Limit, \$\\phi\$",
    title="Interest Rate, \$r\$",
    alpha = 0.2,
    color = default_colors[1],
    ls = :dash,
    label = "Flexible \$r\$"
)
plt_r_Y = plot(plt_r,plt_Y,
    layout = (1,2),
    size = (1200,400)
)
plot!(margin = 6mm)
if fig_save == 1
    savefig(plt_r_Y, fig_dir*"Y_r_ss_counterfactual.pdf")
end

########################################################
# Permanent shock to phi using Sequence Space Jacobian
########################################################
dphi_shock = 0.4
get_IRF_permanent_shock_flexible_r(param_calibrated,dphi_shock)
get_IRF_permanent_shock_rigid_r(param_calibrated,dphi_shock)



########################################################
# Intertemporal MPC
########################################################
r_target = 0.02
phi = 0.0001
a_range = 0.1
J_y = get_iMPC(r_target,phi,a_range)
phi = 1.0
a_range = 5.0;
J_y_loose_phi = get_iMPC(r_target,phi,a_range)

plot_iMPC(param,J_y,J_y_loose_phi)  