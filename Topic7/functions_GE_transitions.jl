function get_IRF_permanent_shock_flexible_r(param_calibrated,dphi_shock)
    @unpack T = param_calibrated
    dphi = 1e-3;
    param_initial_ss = set_parameters(r=r_target,beta=beta_calibrated,transform_phi = true, phi=param_calibrated.phi+dphi)
    r_initial_ss = compute_r(param_initial_ss)
    param_initial_ss = set_parameters(r=r_initial_ss,beta=beta_calibrated,transform_phi = true, phi=param_calibrated.phi + dphi)
    Bellman_result_initial = solve_policy_EGM(param_initial_ss)
    ss_distribution_initial = solve_ss_distribution(param_initial_ss,Bellman_result_initial)
    Cold = Bellman_result_initial.c_pol[:]'*ss_distribution_initial[:]

    #New Steady State
    Bellman_result = solve_policy_EGM(param_calibrated)
    ss_distribution = solve_ss_distribution(param_calibrated,Bellman_result)
    Lambda =  construct_transition_matrix(param_calibrated, Bellman_result.a_pol)
    E = get_expectations(param_calibrated, Lambda, Bellman_result.c_pol)
    Cnew = Bellman_result.c_pol[:]'*ss_distribution[:]

    dr_ss = (r_initial_ss - param_calibrated.r)/dphi
    # Compute the path of consumption from distribution
    Cpath = E[:,1:(T)]'*ss_distribution_initial[:]
    dC_distribution = (Cpath .- Cnew)./dphi


    # Sequence Space Jacobian
    J_rnext, J_y, J_r0,J_phi0, J_phinext = SSJ_wrapper(param_calibrated,Bellman_result, ss_distribution; Outcome = "C")

    dC = (dC_distribution + J_phi0)*dphi_shock
    plot(1:T,dC)
    plot!(1:T,dC_distribution)
    plot!(1:T,J_phi0)
    IRF_r = -J_rnext\dC .- dr_ss*dphi_shock .+ param_calibrated.r

    r_preshock = param_calibrated.r

    Tmax_plt = 30;
    preT = 5;
    plt_idx = 1:(Tmax_plt+preT)
    IRF_r_plot = vcat(ones(preT)*r_preshock,IRF_r)
    t_vec = -preT:(T-1)
    plt_r = plot(t_vec[plt_idx],IRF_r_plot[plt_idx],
        lw = 5,
        label = :none,
        xlabel = "Year",
        title = "Interest Rate, \$r\$"
    )

    IRF_phi = param_calibrated.phi .- dphi_shock*ones(T)
    IRF_phi_plot = vcat(ones(preT)*param_calibrated.phi,IRF_phi)
    plt_phi = plot(t_vec[plt_idx],IRF_phi_plot[plt_idx],
        lw = 5,
        label = :none,
        xlabel = "Year",
        title = "Borrowing Limit, \$\\phi\$"
    )

    IRF_Y_plot = ones(length(t_vec))*param_calibrated.Y
    plt_Y = plot(t_vec[plt_idx],IRF_Y_plot[plt_idx],
        lw = 5,
        label = :none,
        xlabel = "Year",
        title = "Output, \$Y\$"
    )
    ylims!((0.5,1.2))

    plt_all = plot(plt_phi,plt_r,plt_Y,layout = (1,3),size = (1200,400))
    plot!(margin = 5mm)
    if fig_save == 1
        savefig(plt_all,fig_dir*"IRF_flexible_r.pdf")
    end
    return plt_all
end


function get_IRF_permanent_shock_rigid_r(param_calibrated,dphi_shock)
    @unpack T = param_calibrated
    dphi = 1e-3;
    param_initial_ss = set_parameters(r=r_target,beta=beta_calibrated,transform_phi = true, phi=param_calibrated.phi+dphi)
    Y_initial_ss = compute_Y(param_initial_ss)
    param_initial_ss = set_parameters(r = r_target, Y = Y_initial_ss,beta=beta_calibrated,transform_phi = true, phi=param_calibrated.phi + dphi)
    Bellman_result_initial = solve_policy_EGM(param_initial_ss)
    ss_distribution_initial = solve_ss_distribution(param_initial_ss,Bellman_result_initial)
    Cold = Bellman_result_initial.c_pol[:]'*ss_distribution_initial[:]

    #New Steady State
    Bellman_result = solve_policy_EGM(param_calibrated)
    ss_distribution = solve_ss_distribution(param_calibrated,Bellman_result)
    Lambda =  construct_transition_matrix(param_calibrated, Bellman_result.a_pol)
    E = get_expectations(param_calibrated, Lambda, Bellman_result.c_pol)
    Cnew = Bellman_result.c_pol[:]'*ss_distribution[:]

    dY_ss = (Y_initial_ss - param_calibrated.Y)/dphi
    # Compute the path of consumption from distribution
    Cpath = E[:,1:(T)]'*ss_distribution_initial[:]
    dC_distribution = (Cpath .- Cnew)./dphi


    # Sequence Space Jacobian
    J_rnext, J_y, J_r0,J_phi0, J_phinext = SSJ_wrapper(param_calibrated,Bellman_result, ss_distribution; Outcome = "C")

    dC = (dC_distribution + J_phi0)*dphi_shock
    plot(1:T,dC)
    plot!(1:T,dC_distribution)
    plot!(1:T,J_phi0)
    IRF_Y = -J_y\dC .- dY_ss*dphi_shock .+ param_calibrated.Y

    Y_preshock = param_calibrated.Y

    Tmax_plt = 30;
    preT = 5;
    plt_idx = 1:(Tmax_plt+preT)
    IRF_Y_plot = vcat(ones(preT)*Y_preshock,IRF_Y)
    t_vec = -preT:(T-1)
    plt_Y = plot(t_vec[plt_idx],IRF_Y_plot[plt_idx],
        lw = 5,
        label = :none,
        xlabel = "Year",
        title = "Output, \$Y\$"
    )

    IRF_phi = param_calibrated.phi .- dphi_shock*ones(T)
    IRF_phi_plot = vcat(ones(preT)*param_calibrated.phi,IRF_phi)
    plt_phi = plot(t_vec[plt_idx],IRF_phi_plot[plt_idx],
        lw = 5,
        label = :none,
        xlabel = "Year",
        title = "Borrowing Limit, \$\\phi\$"
    )

    IRF_r_plot = ones(length(t_vec))*param_calibrated.r
    plt_r = plot(t_vec[plt_idx],IRF_r_plot[plt_idx],
        lw = 5,
        label = :none,
        xlabel = "Year",
        title = "Interest Rate, \$r\$"
    )
    ylims!((0.0,0.03))


    plt_all = plot(plt_phi,plt_r,plt_Y,layout = (1,3),size = (1200,400))
    plot!(margin = 5mm)
    if fig_save == 1
        savefig(plt_all,fig_dir*"IRF_rigid_r.pdf")
    end
    return plt_all
end

