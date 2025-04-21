using ForwardDiff

function solve_policy_EGM(param, param_changed)
    @unpack Na, Ny,tol = param
    @unpack beta,r,phi,Y = param_changed
    c_pol_old = ones(Na,Ny)
    c_pol_new = 100*ones(Na,Ny)
    a_pol_new = []
    iter = 0;
    while maximum(abs.(c_pol_new .- c_pol_old)) > tol
        c_pol_old = c_pol_new
        c_pol_new, a_pol_new = Euler_iteration_once(param,c_pol_old,param_changed)
        iter += 1
    end
    println("Policy iteration converged in $iter iterations with maximum error: $(maximum(abs.(c_pol_new .- c_pol_old)))")
    return (
        c_pol = c_pol_new,
        a_pol = a_pol_new,
    )
end


function Euler_iteration_once(param,c_pol_old, param_changed)
    @unpack Na, Ny, ag, yg,ytran,amin = param
    @unpack beta,r,phi,Y = param_changed

    uprime_future = uprime_fun(param,c_pol_old)
    a_today_unconstrained = zeros(eltype(uprime_future),Na,Ny) .+ zeros(eltype(r),1)

    for (ia,a_future) in enumerate(ag)
        uprime_future_y =  uprime_future[ia,:]
        Euler_RHS = beta.*(1.0 .+r).*(ytran*uprime_future_y)
        c_today_unconstrained = uprime_inv_fun(param,Euler_RHS)
        a_today_unconstrained[ia,:] =  (c_today_unconstrained .+ a_future .- yg.*Y .+ r.*phi )./(1.0 .+r) .+ zeros(eltype(a_today_unconstrained),Ny)
    end

    a_pol_new = zeros(eltype(a_today_unconstrained),Na,Ny)
    c_pol_new = zeros(eltype(a_today_unconstrained),Na,Ny)
    for (iy,y_today) in enumerate(yg)
        ainterp = LinearInterpolation(a_today_unconstrained[:,iy],ag, extrapolation_bc=Interpolations.Flat())
        a_pol_new[:,iy] = ainterp.(ag)
        a_pol_new[a_pol_new[:,iy] .< amin,iy] .= amin
        c_pol_new[:,iy] = (1 .+r).*ag .+ y_today.*Y .-r.*phi .- a_pol_new[:,iy]
    end
    return c_pol_new, a_pol_new
end 
