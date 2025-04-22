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
