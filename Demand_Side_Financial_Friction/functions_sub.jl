function uprime_fun(param,c)
    @unpack gamma = param        
    return c.^(-gamma)
end


function uprime_inv_fun(param,uprime)
    @unpack gamma = param        
    return uprime.^(-1/gamma)
end

function compute_index_ia_iy(param,ia,iy)
    @unpack Na, Ny = param
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
    else
        share_distance = distances ./ (ag[indices_above] - ag[indices_above-1])
    end
    
    return indices_above, share_distance
end
