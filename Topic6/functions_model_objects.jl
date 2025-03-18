
function solve_for_r(param,omega,lambda)
    @unpack zg = param
   
    integrate_omega_from_above = reverse(cumsum(reverse(omega)))
    capital_market_clearing = lambda .* integrate_omega_from_above .- 1
    r_index = findlast(capital_market_clearing .> 0)
    if isnothing(r_index)
        r_index = 1
    end
    r = zg[r_index]

    return r,r_index
end

function compute_piz(param,r,lambda)
    @unpack zg = param
    piz = max.(zg .- r,0).*(lambda-1) .+ r
    return piz
end

function compute_Z(param,omega,r_index)
    @unpack zg = param
    fraction_above_r = sum(omega[r_index:end])
    int_from_r_zmax = sum(zg[r_index:end].*omega[r_index:end])

    Z = int_from_r_zmax/fraction_above_r
    return Z
end

function compute_Kgrowth(param,omega,piz)
    @unpack zg,beta = param
    Kgrowth = sum(beta.*(1 .+piz).*omega)
    return Kgrowth
end

function compute_omega_next(param,omega,piz,Kgrowth)
    @unpack ztran, beta = param
    saving_from_z = (beta.*(1 .+ piz).*omega)
    omega_next = 1/Kgrowth.* ( ztran'*saving_from_z )
    return omega_next
end
