
function solve_eqm_path(param,omega_init,lambda_path, T)
    @assert length(lambda_path) == T
    @unpack Nz,tol = param

    omega_path = zeros(Nz,T)
    omega_path[:,1] = omega_init;
    r_path = zeros(T)
    r_index_path = zeros(T)
    Z_path = zeros(T)
    Kgrowth_path = zeros(T)

    Tend = []
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
        if t < T
            omega_path[:,t+1] = omega_next;
            omega_diff = maximum(abs.(omega_path[:,t+1] .- omega_path[:,t]));
            println("t = $t, omega_diff = $omega_diff")
            if omega_diff < tol
                println("Steady state reached at t = $t")
                Tend = copy(t);
                break
            else
                Tend = copy(T);
            end
        end        
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
