
function ss_comparative_statics(param, lambda_vec; T = 10000, omega_in = nothing )
    @unpack Nz = param
    Z_vec = zeros(length(lambda_vec))
    omega_vec = zeros(Nz,length(lambda_vec))
    r_vec = zeros(length(lambda_vec))



    for (i_lambda,lambda) in enumerate(lambda_vec)
        if isnothing(omega_in)
            omega_init = ones(Nz)/Nz
        else
            omega_init = copy(omega_in)
        end
        lambda_path = lambda*ones(T);
        eqm_result = solve_eqm_path(param,omega_init,lambda_path, T)
        Z_vec[i_lambda] = eqm_result.Z_path[end]
        omega_vec[:,i_lambda] = eqm_result.omega_path[:,end]
        r_vec[i_lambda] = eqm_result.r_path[end]
    end

    return (
        Z_vec = Z_vec,
        omega_vec = omega_vec,
        r_vec = r_vec
    )
end


function solve_eqm_path(param,omega_init,lambda_path, T)
    @assert length(lambda_path) == T
    @unpack Nz,tol = param

    omega_path = zeros(Nz,T)
    omega_path[:,1] = omega_init;
    r_path = zeros(T)
    r_index_path = zeros(T)
    Z_path = zeros(T)
    Kgrowth_path = zeros(T)

    Tend = T;
    for t = 1:T
        omega = omega_path[:,t];
        lambda = lambda_path[t];
        r, r_index = solve_for_r(param,omega,lambda)
        piz = compute_piz(param,r,lambda)
        Z = compute_Z(param,omega,r_index);
        Kgrowth = compute_Kgrowth(param,omega,piz)
        omega_next = compute_omega_next(param,omega,piz,Kgrowth)

        r_path[t] = r;
        r_index_path[t] = r_index;
        Z_path[t] = Z;
        Kgrowth_path[t] = Kgrowth;
        if t > 1
            Z_diff = abs(Z_path[t] - Z_path[t-1])
            if t % 100 == 0
                println("t = $t, Z_diff = $Z_diff")
            end
            if Z_diff < tol
                println("Steady state reached at t = $t")
                Tend = copy(t);
                break
            end
        end
        if t < T
            omega_path[:,t+1] = omega_next;
        end        
    end

    if Tend < T 
        println("Steady state not reached")
    end

    return (
        omega_path = omega_path[:,1:Tend],
        r_path = r_path[1:Tend],
        r_index_path = r_index_path[1:Tend],
        Z_path = Z_path[1:Tend],
        Kgrowth_path = Kgrowth_path[1:Tend],
        Tend = Tend
    )
end
