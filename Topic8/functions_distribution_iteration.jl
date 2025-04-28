
function solve_ss_distribution(param,Bellman_result)
    @unpack Na, Ny, Nz = param
    @unpack a_pol = Bellman_result
    transition_matrix = construct_transition_matrix(param,a_pol)

    Matrix_to_invert = I - transition_matrix'
    Matrix_to_invert[end,:] = ones(Na*Ny*Nz)

    RHS = zeros(Na*Ny*Nz)
    RHS[end] = 1.0;
    ss_distribution = Matrix_to_invert\RHS
    @assert sum(ss_distribution) ≈ 1.0
    ss_distribution = reshape(ss_distribution,Na,Ny,Nz)
    return ss_distribution
end


function construct_transition_matrix(param, a_pol)
    @unpack Na, Ny, yg, ytran,ag,r,death_prob ,amin,Nz,z_probs = param
    transition_matrix = zeros(eltype(a_pol),Na*Ny*Nz,Na*Ny*Nz)
    for ia = 1:Na
        for iy = 1:Ny
            for iz = 1:Nz
                a_next = a_pol[ia,iy,iz]
                left_grid,right_grid,left_weight,right_weight = find_nearest_grid(ag,a_next)

                a_next_min = amin
                left_grid_min,right_grid_min,left_weight_min,right_weight_min = find_nearest_grid(ag,a_next_min)

                index_ia_iy_iz = compute_index_ia_iy_iz(param,ia,iy,iz)
                for iy_next = 1:Ny
                    for iz_next = 1:Nz 

                        index_ia_iy_iz_next = compute_index_ia_iy_iz(param,left_grid,iy_next,iz_next)
                        transition_matrix[index_ia_iy_iz,index_ia_iy_iz_next] += (1-death_prob)*z_probs[iz_next]*left_weight*ytran[iy,iy_next]

                        index_ia_iy_iz_next = compute_index_ia_iy_iz(param,right_grid,iy_next,iz_next)
                        transition_matrix[index_ia_iy_iz,index_ia_iy_iz_next] += (1-death_prob)*z_probs[iz_next]*right_weight*ytran[iy,iy_next]

                        index_ia_iy_iz_next = compute_index_ia_iy_iz(param,left_grid_min,iy_next,iz_next)
                        transition_matrix[index_ia_iy_iz,index_ia_iy_iz_next] += death_prob*z_probs[iz_next]*left_weight_min*ytran[iy,iy_next]
        
                        index_ia_iy_iz_next = compute_index_ia_iy_iz(param,right_grid_min,iy_next,iz_next)
                        transition_matrix[index_ia_iy_iz,index_ia_iy_iz_next] += death_prob*z_probs[iz_next]*right_weight_min*ytran[iy,iy_next]
                    end
    
                end
            end
        end
    end
    transition_matrix = sparse(transition_matrix)
    @assert sum(transition_matrix,dims=2) ≈ ones(Na*Ny*Nz)
    return transition_matrix
end
