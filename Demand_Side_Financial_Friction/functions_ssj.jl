#SS objects 
function get_ss_objects(param, beta, r)
    Bellman_result = solve_policy_EGM(param, beta, r)
    ss_distribution = solve_ss_distribution(param,Bellman_result)
    @unpack ag,amin = param
    @unpack c_pol, a_pol = Bellman_result
    
    return c_pol, a_pol, ss_distribution

end 

# We need to obtain the expectation vectors (Lemma3)

function get_expectations(param, Lambda, beta, r, y_ss)
    E = zeros(length(y_ss), param.T+1)
    E[ :, 1] = y_ss[:] #assigns the flattened version of y_ss
    for t = 1:param.T
        E[ :, t+1] = Lambda * E[ :, t]
    end
    return E
end 

## Backward iteration step 

#recentering backward objective function

function backward_fn(dx, T, Na, Ny, beta, r, c_ss, ss_distribution)
    K = zeros(eltype(dx), T)
    D = zeros(eltype(dx), length(ss_distribution), T)
    #recentering 
    c_ref = Euler_iteration_once(param,c_ss, beta, r)[1]
    c_recentered = zeros(eltype(dx),Na,Ny) 
    c_pol_t = zeros(eltype(dx),Na,Ny)
    a_pol_t = zeros(eltype(dx),Na,Ny)
    for t in Iterators.reverse(1:T)
        if t == T
            c_pol_t, a_pol_t = Euler_iteration_once(param,c_ss, beta, r+dx) 
        else 
            c_pol_t, a_pol_t = Euler_iteration_once(param,c_recentered, beta, r )
        end 
        Lambda = construct_transition_matrix(param, a_pol_t)
        K[t] = a_pol_t[:]' * ss_distribution[:]
        D[:, t] = Lambda' * ss_distribution[:]

        c_recentered = c_ss .+ (c_pol_t .- c_ref) 
    end 
    return vcat(K, D[:])
end 




function get_curly_Y_D(param, beta, r, c_ss, ss_distribution)
    @unpack T, Na, Ny = param 

    objective_fn(x) = backward_fn(x[1], T, Na, Ny, beta, r, c_ss, ss_distribution)

    derivative = ForwardDiff.jacobian(objective_fn, [r])

    curlyY = reverse(derivative[1:T])  
    curlyD = reverse(reshape(derivative[(T+1):end],(Na*Ny,T)),dims=2)
    return curlyY, curlyD
end 

function get_SSJ(param, curlyY, curlyD, E)
    @unpack T = param
    # fake news matrix 
    F = zeros(T, T)

    F[1, :] = curlyY
    F[2:T, :] = E[:, 1:(T-1)]' * curlyD

    #J is the SSJ matrix 

    J = zeros(T, T)
    J[1, :] = F[1, :]
    J[:, 1] = F[:, 1] 

    for s in 2:T
        for t in 2:T
            J[t, s] = J[t-1, s-1] + F[t,s]
        end 
    end 
    return J, F
end 
 

    