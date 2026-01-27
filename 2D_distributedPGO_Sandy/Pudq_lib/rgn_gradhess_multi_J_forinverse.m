function [rgrad_X, rgnhess_X, F_X, P_X] = rgn_gradhess_multi_J_forinverse(local_graph)
    
    N = length(local_graph.vertices);
    num_intra_edges = size (local_graph.intra_edges, 2);
    num_inter_edges = size (local_graph.inter_edges, 2);

    num_J_triplets = 24*num_intra_edges + 12*num_inter_edges;
    J_i_triplets = zeros(num_J_triplets, 1);
    J_j_triplets = zeros(num_J_triplets, 1);
    J_v_triplets = zeros(num_J_triplets, 1);

    t_J_id = 0;

    num_E_triplets = 3*(num_intra_edges + num_inter_edges);
    E_i_triplets = zeros(num_E_triplets, 1);
    E_j_triplets = ones(num_E_triplets, 1);
    E_v_triplets = zeros(num_E_triplets, 1);

    t_E_id = 0;

    % Intra edges
    for ij=1:num_intra_edges
        edge_i = local_graph.intra_edges(1,ij);
        edge_j = local_graph.intra_edges(2,ij);
        x_i = local_graph.vertices_pudq{edge_i};
        x_j = local_graph.vertices_pudq{edge_j};
        z_ij = local_graph.intra_dp_pudq{ij};

        r_ij = pudq_compose(pudq_inv(z_ij), pudq_mul(pudq_inv(x_i), x_j));

        mu_i    =  z_ij(1)*x_j(1) + z_ij(2)*x_j(2);
        omega_i = -z_ij(2)*x_j(1) + z_ij(1)*x_j(2);
        eta_i   = -z_ij(2)*x_j(1) + z_ij(1)*x_j(2);
        kappa_i = -x_j(1)*z_ij(1) - x_j(2)*z_ij(2);
        alpha_1 = -x_j(1)*z_ij(3) - x_j(2)*z_ij(4) + z_ij(1)*x_j(3)+z_ij(2)*x_j(4);
        beta_1  =  x_j(1)*z_ij(4) - x_j(2)*z_ij(3) - z_ij(2)*x_j(3)+z_ij(1)*x_j(4);
        xi_1    = -x_j(1)*z_ij(1) + z_ij(2)*x_j(2);
        zeta_1  = -z_ij(2)*x_j(1) - z_ij(1)*x_j(2);
        alpha_3 = -x_j(1)*z_ij(4) + x_j(2)*z_ij(3) - x_j(3)*z_ij(2)+z_ij(1)*x_j(4);
        beta_3  = -x_j(1)*z_ij(3) - x_j(2)*z_ij(4) - x_j(3)*z_ij(1)-z_ij(2)*x_j(4);
        
        mu_j    =  x_i(1)*z_ij(1) - x_i(2)*z_ij(2);
        omega_j =  x_i(1)*z_ij(2) + x_i(2)*z_ij(1);
        eta_j   = -z_ij(2)*x_i(1) - z_ij(1)*x_i(2);
        kappa_j =  z_ij(1)*x_i(1) - z_ij(2)*x_i(2);
        alpha_2 = -z_ij(3)*x_i(1) + z_ij(4)*x_i(2) - z_ij(1)*x_i(3) - z_ij(2)*x_i(4);
        beta_2  = -z_ij(4)*x_i(1) - z_ij(3)*x_i(2) + z_ij(2)*x_i(3) - z_ij(1)*x_i(4);
        
        phi = get_phi_atan2(r_ij(2), r_ij(1));
        gamma = sinc1(phi);
        
        eij = [r_ij(2)/gamma; r_ij(3)/gamma; r_ij(4)/gamma];
        f1 = f_1(phi);

        Aij = zeros(3,4);
        Bij = zeros(3,4);

        dphi_dxi0 = eta_i*r_ij(1) - mu_i*r_ij(2);
        dphi_dxi1 = kappa_i*r_ij(1) - omega_i*r_ij(2);

        Aij(1,1) = eta_i/gamma   + r_ij(2)*dphi_dxi0*f1;
        Aij(1,2) = kappa_i/gamma + r_ij(2)*dphi_dxi1*f1;
        Aij(2,1) = alpha_1/gamma + r_ij(3)*dphi_dxi0*f1;
        Aij(2,2) = beta_1/gamma  + r_ij(3)*dphi_dxi1*f1;
        Aij(2,3) = xi_1/gamma;
        Aij(2,4) = zeta_1/gamma;
        Aij(3,1) = alpha_3/gamma + r_ij(4)*dphi_dxi0*f1;
        Aij(3,2) = beta_3/gamma  + r_ij(4)*dphi_dxi1*f1;
        Aij(3,3) = -zeta_1/gamma;
        Aij(3,4) = xi_1/gamma;

        dphi_dxj0 = eta_j*r_ij(1) - mu_j*r_ij(2);
        dphi_dxj1 = kappa_j*r_ij(1) - omega_j*r_ij(2);

        Bij(1,1) = eta_j/gamma + r_ij(2)*dphi_dxj0*f1;
        Bij(1,2) = kappa_j/gamma + r_ij(2)*dphi_dxj1*f1;
        Bij(2,1) = alpha_2/gamma + r_ij(3)*dphi_dxj0*f1;
        Bij(2,2) = beta_2/gamma + r_ij(3)*dphi_dxj1*f1;
        Bij(2,3) = kappa_j/gamma;
        Bij(2,4) = -eta_j/gamma;
        Bij(3,1) = beta_2/gamma + r_ij(4)*dphi_dxj0*f1;
        Bij(3,2) = -alpha_2/gamma + r_ij(4)*dphi_dxj1*f1;
        Bij(3,3) = eta_j/gamma;
        Bij(3,4) = kappa_j/gamma;

        ij_index = 3*(ij-1)+1;

        for i=1:3
            t_E_id = t_E_id+1;
            E_i_triplets(t_E_id) = ij_index+(i-1);
            E_v_triplets(t_E_id) = eij(i);

            for j=1:4
                i_index = 3*(ij-1)+i;
                j_index_A = 4*(edge_i-1)+j;
                j_index_B = 4*(edge_j-1)+j;

                t_J_id = t_J_id+1;
                J_i_triplets(t_J_id) = i_index;
                J_j_triplets(t_J_id) = j_index_A;
                J_v_triplets(t_J_id) = Aij(i,j);

                t_J_id = t_J_id+1;
                J_i_triplets(t_J_id) = i_index;
                J_j_triplets(t_J_id) = j_index_B;
                J_v_triplets(t_J_id) = Bij(i,j);
            end
        end
    end
   
    % Inter edges
    for ij=1:num_inter_edges
        edge_inverse = local_graph.inter_edges(3,ij);
        edge_i = local_graph.inter_edges(1,ij);
        edge_j = local_graph.inter_edges(2,ij);

        ij_index = 3*(num_intra_edges+ij-1)+1;

        if edge_inverse == false
            x_i = local_graph.vertices_pudq{edge_i};
            z_ij = local_graph.inter_dp_pudq{ij};
            x_j = local_graph.lm_vertices_pudq{edge_j};
    
            r_ij = pudq_compose(pudq_inv(z_ij), pudq_mul(pudq_inv(x_i), x_j));
    
            mu_i    =  z_ij(1)*x_j(1) + z_ij(2)*x_j(2);
            omega_i = -z_ij(2)*x_j(1) + z_ij(1)*x_j(2);
            eta_i   = -z_ij(2)*x_j(1) + z_ij(1)*x_j(2);
            kappa_i = -x_j(1)*z_ij(1) - x_j(2)*z_ij(2);
            alpha_1 = -x_j(1)*z_ij(3) - x_j(2)*z_ij(4) + z_ij(1)*x_j(3)+z_ij(2)*x_j(4);
            beta_1  =  x_j(1)*z_ij(4) - x_j(2)*z_ij(3) - z_ij(2)*x_j(3)+z_ij(1)*x_j(4);
            xi_1    = -x_j(1)*z_ij(1) + z_ij(2)*x_j(2);
            zeta_1  = -z_ij(2)*x_j(1) - z_ij(1)*x_j(2);
            alpha_3 = -x_j(1)*z_ij(4) + x_j(2)*z_ij(3) - x_j(3)*z_ij(2)+z_ij(1)*x_j(4);
            beta_3  = -x_j(1)*z_ij(3) - x_j(2)*z_ij(4) - x_j(3)*z_ij(1)-z_ij(2)*x_j(4);
  
            phi = get_phi_atan2(r_ij(2), r_ij(1));
            gamma = sinc1(phi);
            eij = [r_ij(2)/gamma; r_ij(3)/gamma; r_ij(4)/gamma];
            f1 = f_1(phi);
    
            Aij = zeros(3,4);
            dphi_dxi0 = eta_i*r_ij(1) - mu_i*r_ij(2);
            dphi_dxi1 = kappa_i*r_ij(1) - omega_i*r_ij(2);
            Aij(1,1) = eta_i/gamma   + r_ij(2)*dphi_dxi0*f1;
            Aij(1,2) = kappa_i/gamma + r_ij(2)*dphi_dxi1*f1;
            Aij(2,1) = alpha_1/gamma + r_ij(3)*dphi_dxi0*f1;
            Aij(2,2) = beta_1/gamma  + r_ij(3)*dphi_dxi1*f1;
            Aij(2,3) = xi_1/gamma;
            Aij(2,4) = zeta_1/gamma;
            Aij(3,1) = alpha_3/gamma + r_ij(4)*dphi_dxi0*f1;
            Aij(3,2) = beta_3/gamma  + r_ij(4)*dphi_dxi1*f1;
            Aij(3,3) = -zeta_1/gamma;
            Aij(3,4) = xi_1/gamma;
    
            for i=1:3
                t_E_id = t_E_id+1;
                E_i_triplets(t_E_id) = ij_index+(i-1);
                E_v_triplets(t_E_id) = eij(i);

                for j=1:4
                    i_index = 3*(num_intra_edges+ij-1)+i;
                    j_index_A = 4*(edge_i-1)+j;
  
                    t_J_id = t_J_id+1;
                    J_i_triplets(t_J_id) = i_index;
                    J_j_triplets(t_J_id) = j_index_A;
                    J_v_triplets(t_J_id) = Aij(i,j);
                end
            end

        else   
            x_i = local_graph.lm_vertices_pudq{edge_i};
            z_ij = local_graph.inter_dp_pudq{ij};
            x_j = local_graph.vertices_pudq{edge_j};
    
            r_ij = pudq_compose(pudq_inv(z_ij), pudq_mul(pudq_inv(x_i), x_j));
    
            mu_j    =  x_i(1)*z_ij(1) - x_i(2)*z_ij(2);
            omega_j =  x_i(1)*z_ij(2) + x_i(2)*z_ij(1);
            eta_j   = -z_ij(2)*x_i(1) - z_ij(1)*x_i(2);
            kappa_j =  z_ij(1)*x_i(1) - z_ij(2)*x_i(2);
            alpha_2 = -z_ij(3)*x_i(1) + z_ij(4)*x_i(2) - z_ij(1)*x_i(3) - z_ij(2)*x_i(4);
            beta_2  = -z_ij(4)*x_i(1) - z_ij(3)*x_i(2) + z_ij(2)*x_i(3) - z_ij(1)*x_i(4);
            
            phi = get_phi_atan2(r_ij(2), r_ij(1));
            gamma = sinc1(phi);
            eij = [r_ij(2)/gamma; r_ij(3)/gamma; r_ij(4)/gamma];
            f1 = f_1(phi);

            Bij = zeros(3,4);
            dphi_dxj0 = eta_j*r_ij(1) - mu_j*r_ij(2);
            dphi_dxj1 = kappa_j*r_ij(1) - omega_j*r_ij(2);
            Bij(1,1) = eta_j/gamma + r_ij(2)*dphi_dxj0*f1;
            Bij(1,2) = kappa_j/gamma + r_ij(2)*dphi_dxj1*f1;
            Bij(2,1) = alpha_2/gamma + r_ij(3)*dphi_dxj0*f1;
            Bij(2,2) = beta_2/gamma + r_ij(3)*dphi_dxj1*f1;
            Bij(2,3) = kappa_j/gamma;
            Bij(2,4) = -eta_j/gamma;
            Bij(3,1) = beta_2/gamma + r_ij(4)*dphi_dxj0*f1;
            Bij(3,2) = -alpha_2/gamma + r_ij(4)*dphi_dxj1*f1;
            Bij(3,3) = eta_j/gamma;
            Bij(3,4) = kappa_j/gamma;

            for i=1:3
                t_E_id = t_E_id+1;
                E_i_triplets(t_E_id) = ij_index+(i-1);
                E_v_triplets(t_E_id) = eij(i);

                for j=1:4
                    i_index = 3*(num_intra_edges+ij-1)+i;
                    j_index_B = 4*(edge_j-1)+j;
    
                    t_J_id = t_J_id+1;
                    J_i_triplets(t_J_id) = i_index;
                    J_j_triplets(t_J_id) = j_index_B;
                    J_v_triplets(t_J_id) = Bij(i,j);
                end
            end
        end
    end

    E = sparse(E_i_triplets, E_j_triplets, E_v_triplets, num_E_triplets, 1);
    J_mat = sparse(J_i_triplets, J_j_triplets, J_v_triplets, num_E_triplets, 4*N);
    F_X = 0.5*full(E'*local_graph.Omega_ij*E);
    egrad_X = J_mat'*local_graph.Omega_ij*E;
    egnhess_X = J_mat'*local_graph.Omega_ij*J_mat;
    P_X = P_X_N_sparse(local_graph.vertices_pudq);

    rgrad_X = full(P_X*egrad_X);
    rgnhess_X = P_X*egnhess_X*P_X;
    rgnhess_X = tril(rgnhess_X) + tril(rgnhess_X,-1)';
    rgnhess_X = full(rgnhess_X);
end
