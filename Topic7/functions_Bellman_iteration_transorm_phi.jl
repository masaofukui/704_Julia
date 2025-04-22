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
        phi = param.phi
        Y = param.Y
    else
        beta = param_changed.beta
        r = param_changed.r
        phi = param_changed.phi
        Y = param_changed.Y
    end
    while maximum(abs.(c_pol_new .- c_pol_old)) > tol
        c_pol_old = c_pol_new
        c_pol_new, a_pol_new = Euler_iteration_once(param,c_pol_old,beta,r,r,phi,phi,Y)
        iter += 1
    end
    println("Policy iteration converged in $iter iterations with maximum error: $(maximum(abs.(c_pol_new .- c_pol_old)))")
    return (
        c_pol = c_pol_new,
        a_pol = a_pol_new,
    )
end


function Euler_iteration_once(param,c_pol_old, beta,r_t,r_tm1,phi_t,phi_tm1,Y)
    @unpack Na, Ny, ag, yg,ytran,amin = param

    uprime_future = uprime_fun(param,c_pol_old)
    a_today_unconstrained = zeros(eltype(uprime_future),Na,Ny) .+ zeros(eltype(r_t),1) .+ zeros(eltype(r_tm1),1)

    for (ia,a_future) in enumerate(ag)
        uprime_future_y =  uprime_future[ia,:]
        Euler_RHS = beta.*(1.0 .+r_t).*(ytran*uprime_future_y)
        c_today_unconstrained = uprime_inv_fun(param,Euler_RHS)
        a_today_unconstrained[ia,:] =  (c_today_unconstrained .+ a_future .- phi_t .- yg.*Y .+ (1.0 .+ r_tm1).*phi_tm1 )./(1.0 .+r_tm1) .+ zeros(eltype(a_today_unconstrained),Ny)
    end

    a_pol_new = zeros(eltype(a_today_unconstrained),Na,Ny)
    c_pol_new = zeros(eltype(a_today_unconstrained),Na,Ny)
    for (iy,y_today) in enumerate(yg)
        ainterp = LinearInterpolation(a_today_unconstrained[:,iy],ag, extrapolation_bc=Interpolations.Flat())
        a_pol_new[:,iy] = ainterp.(ag)
        a_pol_new[a_pol_new[:,iy] .< amin,iy] .= amin
        c_pol_new[:,iy] = (1 .+ r_tm1).*ag .+ y_today.*Y .- (1.0 .+ r_tm1).*phi_tm1 .- a_pol_new[:,iy] .+ phi_t
    end
    return c_pol_new, a_pol_new
end 
