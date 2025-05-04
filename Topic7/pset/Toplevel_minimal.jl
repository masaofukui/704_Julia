using QuantEcon
using Parameters
using Plots
using Interpolations
using SparseArrays
using LinearAlgebra
using Roots
using Plots.Measures
include("functions_sub.jl")
include("functions_Bellman_iteration.jl")
include("functions_solve_distribution.jl")
fig_save = true
function set_parameters(;r = nothing, 
    beta = nothing, phi = 2.0, 
    w = 1.0,a_range=10.0,omega = 0.0,
    sig = 0.2)
    # income process
    rho = 0.95;
    mu = -(1-rho)*0.5*sig^2
    mc = QuantEcon.rouwenhorst(7,rho,sqrt(1-rho^2)*sig,mu)
    ytran = mc.p

    yg = exp.(mc.state_values)
    yg = yg.*(1+omega).*w
    Ny = length(yg);

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
    )
end


########################################################
# Bellman Iteration
########################################################
beta = 0.97
r = 0.02
w = 1.0;
function compute_aggregate_wealth(r,w,omega,sig)
    A = 0.0
    for pm = [-1,1]
        param = set_parameters(r = r, beta = beta, omega=pm*omega, phi = 0.0,  w = w,a_range=5.0,sig=sig)
        Bellman_result = solve_policy_EGM(param, beta, r)
        ss_distribution = solve_ss_distribution(param,Bellman_result)
        A += sum(ss_distribution.*param.ag)
    end
    return A
end

omega_vec = 0.0:0.05:0.5
A_vec = compute_aggregate_wealth.(r,w,omega_vec,0.2)
plt_omega = plot(omega_vec,A_vec,
    lw = 5,
    xlabel = "Permanent Income Inequality, \$ \\omega \$",
    label = :none,
    ylabel = "Aggregate Wealth",
    title = "Aggregate Wealth"
)
ylims!((0,1))

sig_vec = 0.15:0.05:0.5
A_vec = compute_aggregate_wealth.(r,w,0.0,sig_vec)
plt_sig = plot(sig_vec,A_vec,
    lw = 5,
    xlabel = "Volatility of Income, \$ \\sigma \$",
    label = :none,
    ylabel = "Aggregate Wealth",
    title = "Aggregate Wealth"
)

plt_all = plot(plt_sig,plt_omega,layout=(1,2),size=(1000,400))
plot!(margin=5mm)

if fig_save
    savefig(plt_all,"./Topic7/pset/figure/aggregate_wealth.pdf")
end



function solve_GE(r,omega,sig)
    delta = 0.08;
    alph = 1/3;
    K = ((r + delta)/alph).^(1/(alph-1));
    w = ((r + delta)/alph).^(alph/(alph-1));
    A = compute_aggregate_wealth(r,w,omega,sig)
    return A - K, K, w
end

omega = 0.2;
sig = 0.2;

r_vec_omega = [];
K_vec_omega = [];
w_vec_omega = [];
for omega in omega_vec
    solve_r = x -> solve_GE(x,omega,sig)[1]
    r_ss_eqm = find_zero(solve_r, (0.0, 1/beta-1-1e-5))
    A, K, w = solve_GE(r_ss_eqm,omega,sig)
    push!(K_vec_omega,K)
    push!(w_vec_omega,w)
    push!(r_vec_omega,r_ss_eqm)
end

plt_K = plot(omega_vec,K_vec_omega,
    lw = 5,
    xlabel = "Permanent Income Inequality, \$ \\omega \$",
    label = :none,
    ylabel = "Capital",
    title = "Capital"
)
ylims!(K_vec_omega[1]-0.05,K_vec_omega[end]+0.05)
plt_w = plot(omega_vec,w_vec_omega,
    lw = 5,
    xlabel = "Permanent Income Inequality, \$ \\omega \$",
    label = :none,
    ylabel = "Wage",
    title = "Wage"
)
ylims!(w_vec_omega[1]-0.05,w_vec_omega[end]+0.05)
plt_r = plot(omega_vec,r_vec_omega,
    lw = 5,
    xlabel = "Permanent Income Inequality, \$ \\omega \$",
    label = :none,
    ylabel = "Interest Rate",
    title = "Interest Rate"
)
ylims!(r_vec_omega[1]-0.05,r_vec_omega[end]+0.05)
plt_all_omega = plot(plt_K,plt_w,plt_r,layout=(1,3),size=(1500,400))
plot!(margin=5mm)

if fig_save
    savefig(plt_all_omega,"./Topic7/pset/figure/equilibrium_omega.pdf")
end

r_vec_sig = [];
K_vec_sig = [];
w_vec_sig = [];
for sig in sig_vec
    solve_r = x -> solve_GE(x,omega,sig)[1]
    r_ss_eqm = find_zero(solve_r, (-0.01, 1/beta-1))
    A, K, w = solve_GE(r_ss_eqm,omega,sig)
    push!(K_vec_sig,K)
    push!(w_vec_sig,w)
    push!(r_vec_sig,r_ss_eqm)
end


plt_K = plot(sig_vec,K_vec_sig,
    lw = 5,
    xlabel = "Volatility of Income, \$ \\sigma \$",
    label = :none,
    ylabel = "Capital",
    title = "Capital"
)
plt_w = plot(sig_vec,w_vec_sig,
    lw = 5,
    xlabel = "Volatility of Income, \$ \\sigma \$",
    label = :none,
    ylabel = "Wage",
    title = "Wage"
)
plt_r = plot(sig_vec,r_vec_sig,
    lw = 5,
    xlabel = "Volatility of Income, \$ \\sigma \$",
    label = :none,
    ylabel = "Interest Rate",
    title = "Interest Rate"
)

plt_all_sig = plot(plt_K,plt_w,plt_r,layout=(1,3),size=(1500,400))
plot!(margin=5mm)

if fig_save
    savefig(plt_all_sig,"./Topic7/pset/figure/equilibrium_sig.pdf")
end
