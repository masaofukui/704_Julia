
function plot_power_law(param,dist_a,ag; fig_save = 0, fig_name = "",title_label ="",low_q=0.1,high_q=0.001)
    @unpack yg, Ny = param
    empirical_cdf_a = cumsum(dist_a)
    empirical_cdf_a = empirical_cdf_a./empirical_cdf_a[end]
    ccdf_a = 1.0.-empirical_cdf_a
    idx_a = ccdf_a .<= low_q .&& ccdf_a .>= high_q
    
    # Calculate linear fit
    x = log.(ag[idx_a])
    y = log.(ccdf_a[idx_a])
    A = [ones(length(x)) x] 
    coeffs = A \ y
    slope = coeffs[2]
    intercept = coeffs[1]
    fit_line = intercept .+ slope.*x
    
    # Plot data and fit
    plt_power_law = plot(x, y,
        lw=5,
        label=:none,
        color=default_colors[1]
    )
    plot!(x, fit_line,
        lw=3,
        label=:none,
        linestyle=:dash,
        color=default_colors[2]
    )
    find_x = findfirst(x .>= mean(x))
    range_y = maximum(y)-minimum(y)
    annotate!(x[find_x], y[find_x]+range_y/5,
        text("Slope = $(round(slope,digits=2))", :black, :left, 10, font("Computer Modern")))

    title!(title_label)
    ylabel!("Log CCDF")
    xlabel!("Log value")

    return plt_power_law
end

function plot_power_law_wrapper(param,Bellman_result,ss_distribution; fig_save = 0, fig_name = "")
    @unpack ag,yg,r,psi,eta = param
    @unpack c_pol = Bellman_result

    dist_a = vec(sum(ss_distribution,dims=[2,3]))
    dist_y = vec(sum(ss_distribution,dims=[1,3]))


    plt_a = plot_power_law(param,dist_a,ag,title_label = "Wealth",low_q = 0.1,high_q = 0.001)
    plt_y = plot_power_law(param,dist_y,yg,title_label = "Labor Income",low_q = 0.1,high_q = 0.001)
    plt_r = plot_power_law(param,dist_a,r.*psi.*ag.^(1+eta),title_label = "Capital Income",low_q = 0.1,high_q = 0.001)

    lc = range(minimum(log.(c_pol)),maximum(log.(c_pol)),length=Na)
    c_vec = exp.(lc)

    dist_c = get_distc(param,Bellman_result,ss_distribution,c_vec)
    plt_c = plot_power_law(param,dist_c,c_vec,title_label = "Consumption",low_q = 0.1,high_q = 0.0001)

    plt_all = plot(plt_y,plt_a,plt_r,plt_c,layout = (2,2),size = (1200,800))
    plot!(margin = 4mm)
    if fig_save == 1
        savefig(plt_all, fig_dir*fig_name*".pdf")
    end
    display(plt_all)
end