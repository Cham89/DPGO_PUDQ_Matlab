function [rgrad_X, rgnhess_X, F_X, P_X] = rgn_gradhess_J(G)

    %Set up the problem
    N = length(G.vertices);
    M = length(G.delta_poses);
    
    % Initialize triplets for J construction
    num_J_triplets = 24*M;
    J_i_triplets = zeros(num_J_triplets,1);
    J_j_triplets = zeros(num_J_triplets,1);
    J_v_triplets = zeros(num_J_triplets,1);

    % Initialize J triplet ID
    t_J_id = 0;
    
    % Initialize residual vector
    % E = zeros(3*M,1);
    % E = sparse(3*M,1);
    num_E_triplets = 3*M;
    E_i_triplets = zeros(num_E_triplets,1);
    E_j_triplets = ones(num_E_triplets,1);
    E_v_triplets = zeros(num_E_triplets,1);

    % Initialize E triplet ID
    t_E_id = 0;

    for ij=1:M
        edge_i = G.edges(ij,1);
        edge_j = G.edges(ij,2);

        x_i = G.vertices_pudq{edge_i};
        x_j = G.vertices_pudq{edge_j};
        z_ij = G.delta_poses_pudq{ij};

        %Compute residual and Jacobians
        r_ij = pudq_compose(pudq_inv(z_ij), pudq_mul(pudq_inv(x_i), x_j));
    
        %x_i constants
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
        
        %x_j constants
        mu_j    =  x_i(1)*z_ij(1) - x_i(2)*z_ij(2);
        omega_j =  x_i(1)*z_ij(2) + x_i(2)*z_ij(1);
        eta_j   = -z_ij(2)*x_i(1) - z_ij(1)*x_i(2);
        kappa_j =  z_ij(1)*x_i(1) - z_ij(2)*x_i(2);
        alpha_2 = -z_ij(3)*x_i(1) + z_ij(4)*x_i(2) - z_ij(1)*x_i(3) - z_ij(2)*x_i(4);
        beta_2  = -z_ij(4)*x_i(1) - z_ij(3)*x_i(2) + z_ij(2)*x_i(3) - z_ij(1)*x_i(4);
        
        %atan2 method for computing phi
        phi = get_phi_atan2(r_ij(2), r_ij(1));
        gamma = sinc1(phi);
        
        %Compute residual
        eij = [r_ij(2)/gamma; r_ij(3)/gamma; r_ij(4)/gamma];
        
        %Compute atan2 gradient term
        f1 = f_1(phi);
        
        %Initialize Aij jacobian
        Aij = zeros(3,4);
    
        %Initialize Bij jacobian
        Bij = zeros(3,4);

        %Aij jacobian
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
        
        %Bij jacobian
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

        % J_block_ij = {};
        % J_block_ij.eij = eij;
        % J_block_ij.Aij = Aij;
        % J_block_ij.Bij = Bij;
        % J_blocks{ij} = J_block_ij;

        ij_index = 3*(ij-1)+1;
        % E(ij_index:ij_index+2) = eij;

        % Now, write Aij and Bij jacobians into J
        for i=1:3
            % E(ij_index + (i-1)) = eij(i);
            t_E_id = t_E_id+1;
            E_i_triplets(t_E_id) = ij_index+(i-1);
            E_v_triplets(t_E_id) = eij(i);

            for j=1:4
                i_index = 3*(ij-1)+i;
                j_index_A = 4*(edge_i-1)+j;
                j_index_B = 4*(edge_j-1)+j;
                
                % Set Aij triplets
                t_J_id = t_J_id+1;
                J_i_triplets(t_J_id) = i_index;
                J_j_triplets(t_J_id) = j_index_A;
                J_v_triplets(t_J_id) = Aij(i,j);
                
                % Set Bij triplets
                t_J_id = t_J_id+1;
                J_i_triplets(t_J_id) = i_index;
                J_j_triplets(t_J_id) = j_index_B;
                J_v_triplets(t_J_id) = Bij(i,j);
            end
        end
    end

    E = sparse(E_i_triplets, E_j_triplets, E_v_triplets, 3*M, 1);

    % Set the jacobian J from triplets for ultimate performance
    J_mat = sparse(J_i_triplets, J_j_triplets, J_v_triplets, 3*M, 4*N);

    F_X = 0.5*full(E'*G.Omega*E);
    egrad_X = J_mat'*G.Omega*E;
    egnhess_X = J_mat'*G.Omega*J_mat;
    
    % Compute sparse projection matrix and riemannian grad/hess
    P_X = P_X_N_sparse(G.vertices_pudq);

    % Convert to full at the very end (we've already saved our time)
    rgrad_X = full(P_X*egrad_X);
    rgnhess_X = P_X*egnhess_X*P_X;

    % Make it symmetric!
    rgnhess_X = tril(rgnhess_X) + tril(rgnhess_X,-1)';
    rgnhess_X = full(rgnhess_X);

    % F = 0.5*sum(fij);
    % 
    % %Euclidean Gradient
    % egrad_F = zeros(4*N, 1);
    % 
    % %Approximate Euclidean Hessian (Gauss-Newton)
    % egnhess_F = zeros(4*N, 4*N);
    % 
    % %Now, loop and aggregate data into sparse matrices
    % for ij=1:M
    %     edge = edges(ij,:);
    %     i = edge(1);
    %     j = edge(2);
    % 
    %     %Gradient and hessian block indices
    %     i_index = 4*(i-1)+1;
    %     j_index = 4*(j-1)+1;
    % 
    %     egrad_F(i_index:i_index+3) = egrad_F(i_index:i_index+3) + grad_i(:,ij);
    %     egrad_F(j_index:j_index+3) = egrad_F(j_index:j_index+3) + grad_j(:,ij);
    % 
    %     egnhess_F(i_index:i_index+3,i_index:i_index+3) = egnhess_F(i_index:i_index+3,i_index:i_index+3) + Hii_tilde(:,:,ij);
    %     egnhess_F(i_index:i_index+3,j_index:j_index+3) = Hij_tilde(:,:,ij);
    %     egnhess_F(j_index:j_index+3,i_index:i_index+3) = Hij_tilde(:,:,ij)';
    %     egnhess_F(j_index:j_index+3,j_index:j_index+3) = egnhess_F(j_index:j_index+3,j_index:j_index+3) + Hjj_tilde(:,:,ij);   
    % end
    % 
    % %Euclidean Gradient
    % rgrad_F = zeros(4*N, 1);
    % 
    % %Approximate Euclidean Hessian (Gauss-Newton)
    % rgnhess_F = zeros(4*N, 4*N);
    % 
    % %Compute the TX_M^N projection matrix
    % P_X = P_X_N(X);
    % 
    % for ij=1:M
    %     edge = edges(ij,:);
    %     i = edge(1);
    %     j = edge(2);
    % 
    %     %Gradient and hessian block indices
    %     i_index = 4*(i-1)+1;
    %     j_index = 4*(j-1)+1;
    % 
    %     Pxi = P_X(i_index:i_index+3,i_index:i_index+3);
    %     Pxj = P_X(j_index:j_index+3,j_index:j_index+3);
    % 
    %     rgrad_F(i_index:i_index+3) = remove_zeros(Pxi*egrad_F(i_index:i_index+3));
    %     rgrad_F(j_index:j_index+3) = remove_zeros(Pxj*egrad_F(j_index:j_index+3));
    % 
    %     rgnhess_F(i_index:i_index+3,i_index:i_index+3) = remove_zeros(Pxi*egnhess_F(i_index:i_index+3,i_index:i_index+3)*Pxi);
    %     Hij = remove_zeros(Pxi*Hij_tilde(:,:,ij)*Pxj);
    %     rgnhess_F(i_index:i_index+3,j_index:j_index+3) = Hij;
    %     rgnhess_F(j_index:j_index+3,i_index:i_index+3) = Hij';
    %     rgnhess_F(j_index:j_index+3,j_index:j_index+3) = remove_zeros(Pxj*egnhess_F(j_index:j_index+3,j_index:j_index+3)*Pxj);
    % end
    % 
    % P_X = sparse(P_X);
    % rgrad_F = sparse(rgrad_F);
    % rgnhess_F = sparse(rgnhess_F);
end



