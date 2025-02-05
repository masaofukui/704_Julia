using ForwardDiff

function solve_policy_EGM(param, beta, r)
    @unpack Na, Ny,tol = param
    c_pol_old = ones(Na,Ny)
    c_pol_new = 100*ones(Na,Ny)
    a_pol_new = []
    iter = 0;
    while maximum(abs.(c_pol_new .- c_pol_old)) > tol
        c_pol_old = c_pol_new
        c_pol_new, a_pol_new = Euler_iteration_once(param,c_pol_old,beta,r)
        iter += 1
    end
    println("Policy iteration converged in $iter iterations with maximum error: $(maximum(abs.(c_pol_new .- c_pol_old)))")
    return (
        c_pol = c_pol_new,
        a_pol = a_pol_new,
    )
end


function Euler_iteration_once(param,c_pol_old, beta, r)
    @unpack Na, Ny, ag, yg,ytran,amin = param

    uprime_future = uprime_fun(param,c_pol_old)
    a_today_unconstrained = zeros(eltype(uprime_future),Na,Ny) .+ zeros(eltype(r),1)

    for (ia,a_future) in enumerate(ag)
        uprime_future_y =  uprime_future[ia,:]
        Euler_RHS = beta.*(1.0 .+r).*(ytran*uprime_future_y)
        c_today_unconstrained = uprime_inv_fun(param,Euler_RHS)
        a_today_unconstrained[ia,:] =  (c_today_unconstrained .+ a_future .- yg)./(1.0 .+r) .+ zeros(eltype(a_today_unconstrained),Ny)
    end

    a_pol_new = zeros(eltype(a_today_unconstrained),Na,Ny)
    c_pol_new = zeros(eltype(a_today_unconstrained),Na,Ny)
    for (iy,y_today) in enumerate(yg)
        ainterp = LinearInterpolation(a_today_unconstrained[:,iy],ag, extrapolation_bc=Flat())
        a_pol_new[:,iy] = ainterp.(ag)
        a_pol_new[a_pol_new[:,iy] .< amin,iy] .= amin
        c_pol_new[:,iy] = (1 .+r).*ag .+ y_today .- a_pol_new[:,iy]
    end
    return c_pol_new, a_pol_new
end 

function Euler_iteration_once_derivative(param, c_pol_old, beta, r)
    # Function that takes r as only input for ForwardDiff
    function euler_wrt_r(r_input)
        c_pol, a_pol = Euler_iteration_once(param, c_pol_old, beta, r_input)
        return vcat(c_pol[:], a_pol[:])  # Return flattened arrays
    end
    
    # Compute derivative
    derivative = ForwardDiff.derivative(euler_wrt_r, r)
    
    # Reshape derivatives back to original dimensions
    @unpack Na, Ny = param
    total_elements = Na * Ny
    
    # Split the derivative vector back into c_pol and a_pol derivatives
    dc_dr = reshape(derivative[1:total_elements], (Na, Ny))
    da_dr = reshape(derivative[(total_elements+1):end], (Na, Ny))
    
    return dc_dr, da_dr
end

function test_derivative(x)
    y = x*(1:10)
    a = LinearInterpolation(1:10,y)
    return a(x)
end
