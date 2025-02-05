
function solve_ss_distribution(param,Bellman_result)
    @unpack Na, Ny = param
    @unpack a_pol = Bellman_result
    transition_matrix = construct_transition_matrix(param,a_pol)

    Matrix_to_invert = I - transition_matrix'
    Matrix_to_invert[end,:] = ones(Na*Ny)

    RHS = zeros(Na*Ny)
    RHS[end] = 1.0;
    ss_distribution = Matrix_to_invert\RHS
    @assert sum(ss_distribution) ≈ 1.0
    ss_distribution = reshape(ss_distribution,Na,Ny)
    return ss_distribution
end


function construct_transition_matrix(param, a_pol)
    @unpack Na, Ny, yg, ytran = param
    transition_matrix = zeros(eltype(a_pol),Na*Ny,Na*Ny)

    for ia = 1:Na
        for iy = 1:Ny
            index_ia_iy = compute_index_ia_iy(param,ia,iy)
            a_next = a_pol[ia,iy]
            indices_above, share_distance = compute_closes_grid_and_distance(param,a_next)
            indices_below = max(1,indices_above - 1)
            for iy_next = 1:Ny
                index_ia_iy_next = compute_index_ia_iy(param,indices_above,iy_next)
                transition_matrix[index_ia_iy,index_ia_iy_next] += (1 - share_distance)*ytran[iy,iy_next]

                index_ia_iy_next = compute_index_ia_iy(param,indices_below,iy_next)
                transition_matrix[index_ia_iy,index_ia_iy_next] += share_distance*ytran[iy,iy_next]
            end
        end
    end
    transition_matrix = sparse(transition_matrix)
    @assert sum(transition_matrix,dims=2) ≈ ones(Ny*Na)

    return transition_matrix
end
