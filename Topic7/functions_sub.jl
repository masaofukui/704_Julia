function uprime_fun(param,c)
    @unpack gamma = param        
    return c.^(-gamma)
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


function compute_closes_grid_and_distance(param,a_next)
    @unpack ag = param
    
    # Find the next highest grid point for each policy value
    indices_above = searchsortedfirst(ag,a_next)
    
    # Ensure we don't exceed grid boundaries
    indices_above = min.(indices_above, length(ag))

    # Compute distances to the grid points above
    distances = ag[indices_above] - a_next

    # compute as a share of the distance to the closest grid point
    if indices_above == 1
        share_distance = 1.0;
    elseif indices_above == length(ag) && a_next >= ag[length(ag)]
        share_distance = 0.0;
    else
        share_distance = distances ./ (ag[indices_above] - ag[indices_above-1])
    end
    
    return indices_above, share_distance
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
