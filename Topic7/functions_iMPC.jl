function get_iMPC(r_target,phi,a_range,higher_risk = false)
    param = set_parameters(r=r_target,transform_phi = true, phi=phi,a_range=a_range,higher_risk = true)
    beta_calibrated = calibrate_beta(param,r_target)
    param_calibrated = set_parameters(r=r_target,beta=beta_calibrated,transform_phi = true, phi=phi,a_range=a_range,higher_risk = true)
    
    Bellman_result = solve_policy_EGM(param_calibrated)
    ss_distribution = solve_ss_distribution(param_calibrated,Bellman_result)
    J_rnext, J_y, J_r0,J_phi0, J_phinext = SSJ_wrapper(param_calibrated,Bellman_result, ss_distribution; Outcome = "C")
    return J_y
end



function plot_iMPC(param,J_y,J_y_loose_phi)
    Tplot = 6
    plt_iMPC = plot(0:(Tplot-1),J_y[1:Tplot,1],
        lw=5,
        xlabel="Year",
        title= "Marginal Propensity to Consume",
        label = "Tight Borrowing Limit"
    )




    # Data
    years = 0:5
    fagereng = [0.51, 0.18, 0.11, 0.06, 0.03, 0.04]
    ci_lower = [0.44, 0.13, 0.08, 0.03, 0.01, 0.02]
    ci_upper = [0.58, 0.22, 0.14, 0.09, 0.06, 0.08]

    # Error bars
    yerr_lower = fagereng .- ci_lower
    yerr_upper = ci_upper .- fagereng

    # Create plot
    plot!(
        years, fagereng;
        yerror = (yerr_lower, yerr_upper),
        label = "Data from Fagereng et al. (2021)",
        seriestype = :scatter,
        marker = (:circle, 6),
        line = (:solid, :black),
        color = :black,
        legend = :topright,
        xlabel = "Year"
        #ylabel = "Marginal Propensity to Consume"
    )
    plot!(size = (700, 450))

    if fig_save == 1
        savefig(plt_iMPC, fig_dir*"iMPC_tight_phi.pdf")
    end 


    # Add dotted zero line
    #hline!([0]; linestyle = :dot, color = :gray, label = "")

    # Adjust plot style


    plot!(0:(Tplot-1),J_y_loose_phi[1:Tplot,1],
        lw=5,
        label = "Loose Borrowing Limit",
        ls = :dash
    )
    if fig_save == 1
        savefig(plt_iMPC, fig_dir*"iMPC_loose_phi.pdf")
    end
    plot!(0:(Tplot-1),ones(Tplot).*(1.0 .- 1/(1+r_target)),
        lw=5,
        label = "Representative Agent",
        ls = :dot
    )
    if fig_save == 1
        savefig(plt_iMPC, fig_dir*"iMPC_complete_market.pdf")
    end

    return plt_iMPC
end