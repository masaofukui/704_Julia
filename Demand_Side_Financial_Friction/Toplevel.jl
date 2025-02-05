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
include("functions_Bellman_iteration.jl")
include("functions_distribution_iteration.jl") 
include("functions_ssj.jl")

fig_save = 1

function set_parameters(;beta = nothing, amin = -2)
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
    amin = -2;
    amax = 20;
    ag = range(amin,amax,length=Na);

    # risk aversion for utility function
    gamma = 1.0;
    # step size for numerical differentiation
    dx = 0.01; 
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
        T,
        dx, 
    )
end

param = set_parameters()


beta = 0.98
r = 0.02;
@assert beta*(1+r) < 1
@assert param.amin > - minimum(param.yg)./r
Bellman_result = solve_policy_EGM(param, beta, r)
ss_distribution = solve_ss_distribution(param,Bellman_result)

@unpack ag,amin = param
@unpack c_pol, a_pol = Bellman_result

# Set the default font to Computer Modern
default(fontfamily="Computer Modern")

a_idx = 1:40
p1 = plot(ag[a_idx], c_pol[a_idx,[2,1]], 
    lw=5,
    label=["high income" "low income" ],
    linestyle=[:solid :dash])
title!(L"Consumption Policy Functions, $c(a,y)$")
xlabel!(L"Asset, $a$")

#plot wealth distribution 
p2 = plot(ag, sum(ss_distribution, dims =2), 
    lw=5,
    linestyle=[:solid :dash])
title!("Wealth distribution")

display(p1)
display(p2)

 
#trying out the functions 
c_ss, a_ss, ss_distribution = get_ss_objects(param, beta, r)
Lambda =  construct_transition_matrix(param, a_ss)
E = get_expectations(param, Lambda, beta, r, a_ss)
curlyY, curlyD = get_curly_Y_D(param, beta, r, c_ss, a_ss, ss_distribution)
 J =  get_SSJ(param, curlyY, curlyD, E)[1]
 F =  get_SSJ(param, curlyY, curlyD, E)[2]
