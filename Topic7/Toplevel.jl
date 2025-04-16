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
fig_save = 1
default_colors = palette(:auto);

function set_parameters(;beta = nothing, phi = 2.0, transform_phi = false)
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
    )
end

param = set_parameters(transform_phi = true, phi=2.0)


fig_save = 1

beta = 0.978
r = 0.02;
param_changed = 
    (beta = beta, 
    r = r, 
    phi = param.phi)

fig_name = "beta_"*string(beta)*"_r_"*string(r);
wrapper_PE_policy_plot(param,param_changed,fig_name = fig_name,fig_save = fig_save)
Bellman_result = solve_policy_EGM(param, param_changed)
plt_all = run_simulation(param, Bellman_result; seed_num = 1231,fig_save = fig_save)
ss_distribution = solve_ss_distribution(param,Bellman_result)

plot_distribution(param, Bellman_result, ss_distribution,fig_save = fig_save, fig_name = fig_name)



Bellman_list = Dict{Int8,Any}()
ss_distribution_list = Dict{Int8,Any}()
phi_vec = [2.0, 1.5]
for (i,phi) in enumerate(phi_vec)
    param_changed = 
        (beta = beta, 
        r = r, 
        phi = phi)
    Bellman_result = solve_policy_EGM(param, param_changed)
    Bellman_list[i] = Bellman_result
    ss_distribution = solve_ss_distribution(param,Bellman_result)
    ss_distribution_list[i] = ss_distribution
end

plot_distribution_compare(param, Bellman_list, ss_distribution_list,phi_vec,fig_save = fig_save,fig_name = "phi")


Bellman_list = Dict{Int8,Any}()
ss_distribution_list = Dict{Int8,Any}()
r_vec = [0.02 0.021]
phi_vec = param.phi*ones(length(r_vec))
C_vec = zeros(length(r_vec))
A_vec = zeros(length(r_vec))
excess_demand = zeros(length(r_vec))
for (i,r) in enumerate(r_vec)
    param_changed = 
        (beta = beta, 
        r = r, 
        phi = param.phi)
    Bellman_result = solve_policy_EGM(param, param_changed)
    Bellman_list[i] = Bellman_result
    ss_distribution = solve_ss_distribution(param,Bellman_result)
    ss_distribution_list[i] = ss_distribution
    C_vec[i] = sum(ss_distribution.*Bellman_result.c_pol)
    ag_phi = ag .- param.phi;
    A_vec[i] = sum(ss_distribution.*ag_phi)

    excess_demand[i] = r*A_vec[i] + Y - C_vec[i]
end
plot_distribution_compare(param, Bellman_list, ss_distribution_list,phi_vec,fig_save = fig_save,fig_name = "r")


plot(r_vec,C_vec,
    lw=5,
    label="Consumption",
    xlabel="Interest rate",
    ylabel="Consumption",
    title="Consumption as a function of interest rate")
plot!(r_vec,Y*ones(length(r_vec)))
plot(r_vec,A_vec,
    lw=5,
    label="Asset",
    xlabel="Interest rate",
    ylabel="Asset",
    title="Asset as a function of interest rate")
plot_distribution_compare(param, Bellman_list, ss_distribution_list,phi_vec,fig_save = fig_save)

#trying out the functions 
c_ss, a_ss, ss_distribution = get_ss_objects(param, beta, r)
Lambda =  construct_transition_matrix(param, a_ss)
E = get_expectations(param, Lambda, beta, r, a_ss)
curlyY, curlyD = get_curly_Y_D(param, beta, r, c_ss, a_ss, ss_distribution)
 J =  get_SSJ(param, curlyY, curlyD, E)[1]
 F =  get_SSJ(param, curlyY, curlyD, E)[2]
