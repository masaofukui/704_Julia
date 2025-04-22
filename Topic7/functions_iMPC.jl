function get_iMPC(r_target,phi,a_range)
    param = set_parameters(r=r_target,transform_phi = true, phi=phi,a_range=a_range)
    beta_calibrated = calibrate_beta(param,r_target)
    param_calibrated = set_parameters(r=r_target,beta=beta_calibrated,transform_phi = true, phi=phi,a_range=a_range)
    
    Bellman_result = solve_policy_EGM(param_calibrated)
    ss_distribution = solve_ss_distribution(param_calibrated,Bellman_result)
    J_rnext, J_y, J_r0,J_phi0, J_phinext = SSJ_wrapper(param_calibrated,Bellman_result, ss_distribution; Outcome = "C")
    return J_y
end