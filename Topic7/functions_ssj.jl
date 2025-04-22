

function SSJ_wrapper(param,Bellman_result, ss_distribution; Outcome = "C")
    @unpack c_pol,a_pol = Bellman_result
    Lambda =  construct_transition_matrix(param, a_pol)
    if Outcome == "C"
        E = get_expectations(param, Lambda, c_pol)
    elseif Outcome == "A"
        E = get_expectations(param, Lambda, a_pol)
    end
    curly_Y_D = get_curly_Y_D(param, Bellman_result, ss_distribution; Outcome = Outcome)


    @unpack curlyY_r, curlyY_y, curlyD_r, curlyD_y, curlyY_phi, curlyD_phi = curly_Y_D

    J_r     =  get_SSJ(param, curlyY_r, curlyD_r, E)[1]
    J_y     =  get_SSJ(param, curlyY_y, curlyD_y, E)[1]
    J_phi   =  get_SSJ(param, curlyY_phi, curlyD_phi, E)[1]
    # Separately construct J_r0 and J_r
    J_r0 = J_r[:,1]
    J_rnext = hcat(J_r[:,2:end], [J_r[1,end];J_r[1:(end-1),end]])
    J_phi0 = J_phi[:,1]
    J_phinext = hcat(J_phi[:,2:end], [J_phi[1,end];J_phi[1:(end-1),end]])
    return (
        J_rnext = J_rnext,
        J_y = J_y,
        J_r0 = J_r0,
        J_phi0 = J_phi0,
        J_phinext = J_phinext
    )
end



# We need to obtain the expectation vectors (Lemma3)

function get_expectations(param, Lambda, y_ss)
    E = zeros(length(y_ss), param.T+1)
    E[ :, 1] = y_ss[:] #assigns the flattened version of y_ss
    for t = 1:param.T
        E[ :, t+1] = Lambda * E[ :, t]
    end
    return E
end 

## Backward iteration step 

#recentering backward objective function
function backward_fn(dx,  param, Bellman_result, ss_distribution; Outcome = "C")
    @unpack T, Na, Ny = param
    @unpack beta, r, Y, phi = param
    c_pol_ss = Bellman_result.c_pol
    dx_r = dx[1]
    dx_y = dx[2]
    dx_phi = dx[3]
    C = zeros(eltype(dx), T)
    A = zeros(eltype(dx), T)
    D = zeros(eltype(dx), length(ss_distribution), T)
    #recentering 
    c_pol_ref = Euler_iteration_once(param, c_pol_ss, beta, r, r,phi,phi, Y)[1]
    c_pol_recentered = zeros(eltype(dx),Na,Ny)  
    c_pol_t = zeros(eltype(dx),Na,Ny)
    a_pol_t = zeros(eltype(dx),Na,Ny)
   
    for t in reverse(1:T)
        if t == T
            phi_t = phi + dx_phi
            c_pol_t, a_pol_t = Euler_iteration_once(param,c_pol_ss, beta, r, r+dx_r, phi, phi+dx_phi, Y+dx_y)
        elseif  t == T-1
            phi_t = phi;
            c_pol_t, a_pol_t = Euler_iteration_once(param,c_pol_recentered, beta, r+dx_r, r, phi+dx_phi, phi ,Y)  
        else 
            phi_t = phi;
            c_pol_t, a_pol_t = Euler_iteration_once(param,c_pol_recentered, beta, r, r, phi, phi,Y)
        end 

        Lambda = construct_transition_matrix(param, a_pol_t)
        C[t] = c_pol_t[:]' * ss_distribution[:]
        A[t] = (a_pol_t[:] .- phi_t)' * ss_distribution[:]
        D[:, t] = Lambda' * ss_distribution[:]

        c_pol_recentered = c_pol_ss .+ (c_pol_t .- c_pol_ref) 
    end 
    if Outcome == "C"
        return vcat(C, D[:])
    elseif Outcome == "A"
        return vcat(A, D[:])
    end
end 


function get_curly_Y_D(param, Bellman_result, ss_distribution; Outcome = "C")
    @unpack T, Na, Ny = param 

    objective_fn(x) = backward_fn(x, param, Bellman_result, ss_distribution; Outcome = Outcome)

    #use autimatic differentiation option 
    jacob = ForwardDiff.jacobian(objective_fn, zeros(3))

    #we need to reverse the outputs, since we performed backward iteration

    curlyY_r = reverse(jacob[1:T,1])  
    curlyY_y = reverse(jacob[1:T,2])  
    curlyY_phi = reverse(jacob[1:T,3])
    curlyD_r = reverse(reshape(jacob[(T+1):end,1],(Na*Ny,T)),dims=2)
    curlyD_y = reverse(reshape(jacob[(T+1):end,2],(Na*Ny,T)),dims=2)
    curlyD_phi = reverse(reshape(jacob[(T+1):end,3],(Na*Ny,T)),dims=2)

    return (
        curlyY_r = curlyY_r,
        curlyY_y = curlyY_y,
        curlyY_phi = curlyY_phi,
        curlyD_r = curlyD_r,
        curlyD_y = curlyD_y,
        curlyD_phi = curlyD_phi
    )
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
 

function Jacobina_sanity_check(param,J_Y_collect,J_r_collect,rss;svec = [1,10])
    @unpack T = param

    cumJy = zeros(length(svec))
    cumJr = zeros(length(svec))
    for (i,s) in enumerate(svec)   
        for t = 1:T
            cumJy[i] += J_Y_collect[t,s]./(1+rss)^(t-s)
            cumJr[i] += J_r_collect[t,s]./(1+rss)^(t-s)
        end
    end
    
    @assert maximum(abs.(cumJy .- 1.0)) < 1e-8
    @assert maximum(abs.(cumJr)) < 1e-8
    return cumJy, cumJr
end
