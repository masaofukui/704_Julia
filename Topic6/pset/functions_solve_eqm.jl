
function ss_comparative_statics(param, rho_vec; T = 10000, omega_in = nothing )
    @unpack Nz = param
    Z_vec = zeros(length(rho_vec))
    omega_vec = zeros(Nz,length(rho_vec))
    r_vec = zeros(length(rho_vec))
    Y_vec = zeros(length(rho_vec))

    for (i_rho,rho) in enumerate(rho_vec)
        if isnothing(omega_in)
            omega_init = zeros(Nz)
            omega_init[1] = 1.0
        else
            omega_init = copy(omega_in)
        end
        param = set_parameters(rho_z = rho)
        eqm_result = solve_eqm_path(param,omega_init, T)
        Z_vec[i_rho] = eqm_result.Z_path[end]
        omega_vec[:,i_rho] = eqm_result.omega_path[:,end]
        r_vec[i_rho] = eqm_result.r_path[end]
        Y_vec[i_rho] = eqm_result.Y_path[end]
    end

    return (
        Z_vec = Z_vec,
        omega_vec = omega_vec,
        r_vec = r_vec,
        Y_vec = Y_vec
    )
end


function solve_eqm_path(param,omega_init, T; steady_state = true)
    @unpack Nz,tol,alph,L,delta,beta,zg,lambda = param

    omega_path = zeros(Nz,T)
    omega_path[:,1] = omega_init;
    r_path = zeros(T)
    underz_index_path = zeros(T)
    Z_path = zeros(T)
    K_path = ones(T)
    kz_path = zeros(Nz,T)
    w_path = zeros(T)
    pid_path = zeros(T)
    Y_path = zeros(T)
    Kgrowth_path = zeros(T)
    underz_path = zeros(T)
    Tend = T;
    for t = 1:T
        K = K_path[t];
        omega = omega_path[:,t];
        underz, underz_index = solve_for_underz(param,omega,lambda)
        Z = compute_Z(param,omega,underz_index);
        w = (1-alph)*Z.*K.^(1-alph).*L.^(-alph);
        pid = alph.*(1-alph).^((1-alph)/alph).*w.^(-(1-alph)/alph);
        r = pid.*underz;
        Y = Z.*K.^alph.*L.^(1-alph);
        K_next = beta.*alph.*Y + beta.*(1-delta).*K;
        Kgrowth = K_next./K;
        sz = beta.*(max.(pid.*zg .- r,0) .+ (1+r-delta)) ;
        omega_next = compute_omega_next(param,omega,sz,Kgrowth)

        r_path[t] = r;
        underz_index_path[t] = underz_index;
        underz_path[t] = underz;
        Y_path[t] = Y;
        Z_path[t] = Z;
        pid_path[t] = pid;
        w_path[t] = w;
        if t > 2
            Z_diff = abs(Z_path[t] - Z_path[t-1])
            K_diff = abs(K_path[t] - K_path[t-1])
            if t % 100 == 0
                println("t = $t, Z_diff = $Z_diff")
            end
            if Z_diff < tol && K_diff < tol
                println("Steady state reached at t = $t")
                Tend = copy(t);
                break
            end
        end
        if t < T
            K_path[t+1] = K_next;
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
            Z_path = Z_path[1:Tend],
            w_path = w_path[1:Tend],
            pid_path = pid_path[1:Tend],
            Y_path = Y_path[1:Tend],
            Kgrowth_path = Kgrowth_path[1:Tend],
            underz_path = underz_path[1:Tend],
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
            Z_path = Z_path,
            Y_path = Y_path,
            K_path = K_path,
            Kgrowth_path = Kgrowth_path,
            Tend = Tend
        )
    end
end
