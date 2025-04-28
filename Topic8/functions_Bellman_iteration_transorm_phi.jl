using ForwardDiff

function solve_policy_EGM(param; param_changed = nothing)
    @unpack Na, Ny,tol = param
    c_pol_old = ones(Na,Ny)
    c_pol_new = 100*ones(Na,Ny)
    a_pol_new = []
    iter = 0;
    if isnothing(param_changed)
        beta = param.beta
        r = param.r
        Y = param.Y
    else
        beta = param_changed.beta
        r = param_changed.r
        Y = param_changed.Y
    end
    while maximum(abs.(c_pol_new .- c_pol_old)) > tol
        c_pol_old = c_pol_new
        c_pol_new, a_pol_new = Euler_iteration_once(param,c_pol_old,beta,r,Y)
        iter += 1
        println("Policy iteration $iter: $(maximum(abs.(c_pol_new .- c_pol_old)))")
    end
    println("Policy iteration converged in $iter iterations with maximum error: $(maximum(abs.(c_pol_new .- c_pol_old)))")
    return (
        c_pol = c_pol_new,
        a_pol = a_pol_new,
    )
end


function Euler_iteration_once(param,c_pol_old, beta,r,Y)
    @unpack Na, Ny, ag, yg,ytran,amin,zg,z_probs,amax = param

    x_today_unconstrained = zeros(eltype(c_pol_old),Na,Ny) .+ zeros(eltype(r),1) .+ zeros(eltype(r),1)

    ag_temp = copy(ag)./(1.0 .+ r)
    for (ia,a_future) in enumerate(ag_temp)
        Euler_RHS = zeros(eltype(c_pol_old),Ny)
        for (iy,y) in enumerate(yg)
            c_interp = LinearInterpolation(ag,c_pol_old[:,iy], extrapolation_bc=Interpolations.Flat())
            for (iz,z) in enumerate(zg)
                future_x = (1.0 .+r.*z)*a_future;
                c_future = c_interp(future_x)
                Euler_RHS += z_probs[iz].*beta.*(1.0 .+r.*z).*(ytran[:,iy]*uprime_fun(param,c_future))
            end
        end
        c_today_unconstrained = uprime_inv_fun(param,Euler_RHS)
        x_today_unconstrained[ia,:] =  (c_today_unconstrained .+ a_future  .- yg.*Y )
    end

    a_pol_new = zeros(eltype(x_today_unconstrained),Na,Ny)
    c_pol_new = zeros(eltype(x_today_unconstrained),Na,Ny)
    for (iy,y_today) in enumerate(yg)
        ainterp = LinearInterpolation(x_today_unconstrained[:,iy],ag_temp, extrapolation_bc=Interpolations.Linear())
        a_pol_new[:,iy] = ainterp.(ag)
        a_pol_new[a_pol_new[:,iy] .< amin,iy] .= amin
        c_pol_new[:,iy] = ag .+ y_today.*Y .- a_pol_new[:,iy]
    end
    return c_pol_new, a_pol_new
end 
