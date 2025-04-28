
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
    @unpack Na, Ny, yg, ytran,ag,r,death_prob ,amin = param
    transition_matrix = zeros(eltype(a_pol),Na*Ny,Na*Ny)
    for ia = 1:Na
        for iy = 1:Ny
            index_ia_iy = compute_index_ia_iy(param,ia,iy)
            for iy_next = 1:Ny
                for iz = 1:Nz 
                    a_next = (1+r*zg[iz])*a_pol[ia,iy]
                    left_grid,right_grid,left_weight,right_weight = find_nearest_grid(ag,a_next)

                    index_ia_iy_next = compute_index_ia_iy(param,left_grid,iy_next)
                    transition_matrix[index_ia_iy,index_ia_iy_next] += (1-death_prob)*z_probs[iz]*left_weight*ytran[iy,iy_next]

                    index_ia_iy_next = compute_index_ia_iy(param,right_grid,iy_next)
                    transition_matrix[index_ia_iy,index_ia_iy_next] += (1-death_prob)*z_probs[iz]*right_weight*ytran[iy,iy_next]
                end
                a_next = amin
                left_grid,right_grid,left_weight,right_weight = find_nearest_grid(ag,a_next)

                index_ia_iy_next = compute_index_ia_iy(param,left_grid,iy_next)
                transition_matrix[index_ia_iy,index_ia_iy_next] += death_prob*left_weight*ytran[iy,iy_next]

                index_ia_iy_next = compute_index_ia_iy(param,right_grid,iy_next)
                transition_matrix[index_ia_iy,index_ia_iy_next] += death_prob*right_weight*ytran[iy,iy_next]

            end
        end
    end
    transition_matrix = sparse(transition_matrix)
    @assert sum(transition_matrix,dims=2) ≈ ones(Ny*Na)
    return transition_matrix
end
