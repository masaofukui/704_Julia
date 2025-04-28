

function solve_policy_EGM_z(param, beta, r)
    @unpack Na, Ny,tol,Nz = param
    c_pol_old = ones(Na,Ny,Nz)
    c_pol_new = 100*ones(Na,Ny,Nz)
    a_pol_new = []
    iter = 0;
    while maximum(abs.(c_pol_new .- c_pol_old)) > tol
        c_pol_old = c_pol_new
        c_pol_new, a_pol_new = Euler_iteration_once_z(param,c_pol_old,beta,r)
        iter += 1
    end
    println("Policy iteration converged in $iter iterations with maximum error: $(maximum(abs.(c_pol_new .- c_pol_old)))")
    return (
        c_pol = c_pol_new,
        a_pol = a_pol_new,
    )
end


function Euler_iteration_once_z(param,c_pol_old, beta, r)
    @unpack Na, Ny, ag, yg,ytran,amin,zg,z_probs,Nz = param

    uprime_future = uprime_fun(param,c_pol_old)
    a_today_unconstrained = zeros(eltype(uprime_future),Na,Ny,Nz) .+ zeros(eltype(r),1)

    for (ia,a_future) in enumerate(ag)
        Euler_RHS = zeros(eltype(uprime_future),Ny)
        for (iz,z) in enumerate(zg) 
            uprime_future_y =  uprime_future[ia,:,iz]
            Euler_RHS += z_probs[iz].*beta.*(1.0 .+r.*z).*(ytran*uprime_future_y)
        end
        c_today_unconstrained = uprime_inv_fun(param,Euler_RHS)
        for (iz_today,z_today) in enumerate(zg)
            a_today_unconstrained[ia,:,iz_today] =  (c_today_unconstrained .+ a_future .- yg)./(1.0 .+r.*z_today) 
        end
    end

    a_pol_new = zeros(eltype(a_today_unconstrained),Na,Ny,Nz)
    c_pol_new = zeros(eltype(a_today_unconstrained),Na,Ny,Nz)
    for (iy,y_today) in enumerate(yg)
        for (iz_today,z_today) in enumerate(zg)
            ainterp = LinearInterpolation(a_today_unconstrained[:,iy,iz_today],ag, extrapolation_bc=Interpolations.Flat())
            a_pol_new[:,iy,iz_today] = ainterp.(ag)
            a_pol_new[a_pol_new[:,iy,iz_today] .< amin,iy,iz_today] .= amin
            c_pol_new[:,iy,iz_today] = (1 .+r.*z_today).*ag .+ y_today .- a_pol_new[:,iy,iz_today]
        end
    end
    return c_pol_new, a_pol_new
end 
