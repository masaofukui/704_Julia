function uprime_fun(param,c)
    @unpack gamma = param        
    return c.^(-gamma)
end

function discretize_assets(amin, amax, n_a)
    # find maximum ubar of uniform grid corresponding to desired maximum amax of asset grid
    ubar = log(1 + log(1 + amax - amin))
    
    # make uniform grid
    u_grid = range(0, ubar, length=n_a)
    
    # double-exponentiate uniform grid and add amin to get grid from amin to amax
    return amin .+ exp.(exp.(u_grid) .- 1) .- 1
end

function uprime_inv_fun(param,uprime)
    @unpack gamma = param        
    return uprime.^(- 1 ./gamma)
end

function compute_index_ia_iy(param,ia,iy)
    @unpack Na, Ny = param
    @assert 1 <= ia <= Na
    @assert 1 <= iy <= Ny
    return ia + (iy-1)*Na
end

function compute_index_ia_iy_iz(param,ia,iy,iz)
    @unpack Na, Ny, Nz = param
    @assert 1 <= ia <= Na
    @assert 1 <= iy <= Ny
    @assert 1 <= iz <= Nz
    return ia + (iy-1)*Na + (iz-1)*Na*Ny
end

function find_nearest_grid(agrid,a)
    a_min = minimum(agrid)
    a_max = maximum(agrid)
    na = length(agrid)
    if a <= a_min
        left_grid = 1
        right_grid = 1
        left_weight = 1
        right_weight = 0
    elseif a >= a_max
        left_grid = na
        right_grid = na
        left_weight = 0
        right_weight = 1
    else
        ia = findfirst(agrid .>= a)
        @assert ia > 1
        if ia > 1
            diffa = agrid[ia] - agrid[ia-1];
            left_weight = (agrid[ia] - a)/diffa;
            right_weight = (a - agrid[ia-1])/diffa
            left_grid = ia-1;
            right_grid = ia;
        end
    end
    
    return left_grid,right_grid,left_weight,right_weight
end
function get_distc(param,Bellman_result,ss_distribution,c_vec)
    @unpack c_pol = Bellman_result
    dist_c = zeros(length(c_vec))
    for (ic,c) in enumerate(c_vec)
        idx_c = c_pol .< c
        dist_c[ic] = sum(ss_distribution[idx_c])
    end
    dist_c = dist_c./dist_c[end]
    return dist_c
end