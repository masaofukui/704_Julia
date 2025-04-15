
function ss_comparative_statics(param, lambda_vec; T = 10000, omega_in = nothing )
    @unpack Nz = param
    Z_vec = zeros(length(lambda_vec))
    omega_vec = zeros(Nz,length(lambda_vec))
    r_vec = zeros(length(lambda_vec))
    kz_vec = zeros(Nz,length(lambda_vec))


    for (i_lambda,lambda) in enumerate(lambda_vec)
        if isnothing(omega_in)
            omega_init = zeros(Nz)
            omega_init[1] = 1.0
        else
            omega_init = copy(omega_in)
        end
        lambda_path = lambda*ones(T);
        eqm_result = solve_eqm_path(param,omega_init,lambda_path, T)
        Z_vec[i_lambda] = eqm_result.Z_path[end]
        omega_vec[:,i_lambda] = eqm_result.omega_path[:,end]
        r_vec[i_lambda] = eqm_result.r_path[end]
        kz_vec[:,i_lambda] = eqm_result.kz_path[:,end]
    end

    return (
        Z_vec = Z_vec,
        omega_vec = omega_vec,
        r_vec = r_vec,
        kz_vec = kz_vec
    )
end


function solve_eqm_path(param,omega_init,lambda_path, T; steady_state = true)
    @assert length(lambda_path) == T
    @unpack Nz,tol = param

    omega_path = zeros(Nz,T)
    omega_path[:,1] = omega_init;
    r_path = zeros(T)
    r_index_path = zeros(T)
    Z_path = zeros(T)
    Kgrowth_path = zeros(T)
    kz_path = zeros(Nz,T)
    Tend = T;
    for t = 1:T
        omega = omega_path[:,t];
        lambda = lambda_path[t];
        r, r_index = solve_for_r(param,omega,lambda)
        piz = compute_piz(param,r,lambda)
        kz = compute_kz(param,r,lambda,omega)
        Z = compute_Z(param,omega,r_index);
        Kgrowth = compute_Kgrowth(param,omega,piz)
        omega_next = compute_omega_next(param,omega,piz,Kgrowth)

        r_path[t] = r;
        r_index_path[t] = r_index;
        Z_path[t] = Z;
        Kgrowth_path[t] = Kgrowth;
        kz_path[:,t] = kz;
        if t > 2
            Z_diff = abs(Z_path[t] - Z_path[t-1])
            if t % 100 == 0
                println("t = $t, Z_diff = $Z_diff")
            end
            if Z_diff < tol || (Z_path[t] < Z_path[t-1] && Z_path[t-1] > Z_path[t-2])
                println("Steady state reached at t = $t")
                Tend = copy(t);
                break
            end
        end
        if t < T
            omega_path[:,t+1] = omega_next;
        end        
    end

    if Tend == T 
        println("Steady state not reached")
    end

    if steady_state
        return (
            omega_path = omega_path[:,1:Tend],
            kz_path = kz_path[:,1:Tend],
            r_path = r_path[1:Tend],
            r_index_path = r_index_path[1:Tend],
            Z_path = Z_path[1:Tend],
            Kgrowth_path = Kgrowth_path[1:Tend],
            Tend = Tend
        )
    else
        for t in (Tend+1):T
            omega_path[:,t] = omega_path[:,Tend]
            kz_path[:,t] = kz_path[:,Tend]
            r_path[t] = r_path[Tend]
            r_index_path[t] = r_index_path[Tend]
            Z_path[t] = Z_path[Tend]
            Kgrowth_path[t] = Kgrowth_path[Tend]
        end
        return (
            omega_path = omega_path,
            kz_path = kz_path,
            r_path = r_path,
            r_index_path = r_index_path,
            Z_path = Z_path,
            Kgrowth_path = Kgrowth_path,
            Tend = Tend
        )
    end
end
