using Random
using Interpolations
using Plots
using LaTeXStrings
using Plots.PlotMeasures
function run_simulation(param, Bellman_result; seed_num = 1234, fig_save = 0)
    @unpack c_pol, a_pol = Bellman_result
    @unpack ytran, ag, yg, Na, Ny = param

    T = 100
    default(fontfamily="Computer Modern")
    # Set random seed for reproducibility
    Random.seed!(seed_num)
    
    # Initialize arrays to store simulation results
    c_sim = zeros(T)
    a_sim = zeros(T+1) 
    y_sim = zeros(T)
    
    # Set initial conditions
    a_sim[1] = ag[2] # Start at middle of asset grid
    y_state = 1 # Start in low income state
    
    # Create interpolation objects for policies
    c_interp = [LinearInterpolation(ag, c_pol[:,iy], extrapolation_bc=Flat()) for iy in 1:Ny]
    a_interp = [LinearInterpolation(ag, a_pol[:,iy], extrapolation_bc=Flat()) for iy in 1:Ny]
    
    # Run simulation
    for t in 1:T
        # Get current income
        y_sim[t] = yg[y_state]
        
        # Interpolate consumption and next period assets
        c_sim[t] = c_interp[y_state](a_sim[t])
        a_sim[t+1] = a_interp[y_state](a_sim[t])
        
        # Draw next period's income state
        y_state = rand() < ytran[y_state,1] ? 1 : 2
    end

    plt_c = plot(1:T, c_sim, lw=3, label=:none)
    title!("Consumption", fontsize=24)
    xlabel!(L"Time, $t$")

    plt_y = plot(1:T, y_sim, lw=3, label=:none,color=:green)
    title!("Income", fontsize=24)
    xlabel!(L"Time, $t$")

    plt_a = plot(1:T, a_sim[1:T], lw=3, label=:none,color=:red)
    title!("Assets", fontsize=24)
    xlabel!(L"Time, $t$")

    plt_all = plot(plt_y,plt_c, plt_a, layout=(3,1), size=(1000, 700))
    plot!(margin = 3mm)
    if fig_save == 1
        savefig(plt_all, "./figure/simulation_results"*string(seed_num)*".pdf")
    end
    return plt_all
end