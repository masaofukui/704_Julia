using Roots

function compute_asset_demand_vec(param,var_input_vec; var_change = "r")
    A_vec = zeros(length(var_input_vec))
    Bellman_list = Dict{Int8,Any}()
    ss_distribution_list = Dict{Int8,Any}()
    for (i,var_input) in enumerate(var_input_vec)
        if var_change == "r"
            param_changed = 
                (beta = param.beta, 
                r = var_input, 
                phi = param.phi,
                Y = 1.0)
        elseif var_change == "phi"
            param_changed = 
                (beta = param.beta, 
                r = param.r, 
                phi = var_input,
                Y = 1.0)
        elseif var_change == "Y"
            param_changed = 
                (beta = param.beta, 
                r = param.r, 
                phi = param.phi,
                Y = var_input)
        end
        Bellman_result = solve_policy_EGM(param, param_changed)
        Bellman_list[i] = Bellman_result
        ss_distribution = solve_ss_distribution(param,Bellman_result)
        ss_distribution_list[i] = ss_distribution
        ag_phi = param.ag .- param.phi;
        A_vec[i] = sum(ss_distribution.*ag_phi)
    end
    return A_vec, Bellman_list, ss_distribution_list
end

function compute_asset_demand(param,beta,r)
    param_changed = 
        (beta = beta, 
        r = r,
        phi = param.phi,
        Y = param.Y)
    Bellman_result = solve_policy_EGM(param, param_changed)
    ss_distribution = solve_ss_distribution(param,Bellman_result)
    ag_phi = param.ag .- param.phi;
    A = sum(ss_distribution.*ag_phi)
    return A
end

function calibrate_beta(param,rtarget)
    solve_beta = x -> compute_asset_demand(param,x,rtarget) - 0.0
    beta = find_zero(solve_beta, (0.7, 0.99))
    return beta
end

function compute_r(param)
    solve_r = x -> compute_asset_demand(param,param.beta,x) - 0.0
    r = find_zero(solve_r, (-0.1, 1/param.beta-1-1e-5))
    return r
end