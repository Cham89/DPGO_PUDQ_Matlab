function [RPE, linear, angular] = RPE_G_eucl(graph)
    RPE = 0.0;
    linear_error = 0;
    angular_error = 0;
    
    sz_edges = size(graph.delta_poses);
    num_edges = sz_edges(2);
    for ij=1:num_edges
        z_ij_true = graph.delta_poses_true{ij};

        x_i = graph.vertices{graph.edges(ij,1)};
        x_j = graph.vertices{graph.edges(ij,2)};

        R_i = R_theta(x_i(3));

        dt = R_i' * (x_j(1:2) - x_i(1:2));
        dth = wrap_theta(x_j(3) - x_i(3));
        z_ij_hat = [dt', dth']';

        eij = zeros(3,1);
        eij(1:2) = z_ij_true(1:2)-z_ij_hat(1:2);
        eij(3) = theta_diff(z_ij_true(3), z_ij_hat(3));

        linear_error = linear_error + eij(1:2)'*eij(1:2);
        angular_error = angular_error + eij(3)^2;

        RPE = RPE + eij' * eij;
    end
    
    linear = sqrt(linear_error/num_edges);
    angular = sqrt(angular_error/num_edges);
    
    RPE = sqrt(RPE/num_edges);
end
