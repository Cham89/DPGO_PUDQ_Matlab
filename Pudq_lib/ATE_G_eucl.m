%Absolute trajectory error
function [ATE, linear, angular] = ATE_G_eucl(graph)
    ATE = 0.0;
    linear_error = 0.0;
    angular_error = 0.0;
    
    num_vertices = length(graph.vertices);
    for i=1:num_vertices
        x_i_true = graph.vertices_true{i};
        x_i_hat = graph.vertices{i};
        
        eij = zeros(3,1);
        eij(1:2) = x_i_true(1:2)-x_i_hat(1:2);
        eij(3) = theta_diff(x_i_true(3), x_i_hat(3));

        % t_err = sqrt(eij(1:2)'*eij(1:2))

        linear_error = linear_error + eij(1:2)'*eij(1:2);
        angular_error = angular_error + eij(3)^2;

        ATE = ATE + eij' * eij;
    end
    
    linear = sqrt(linear_error/num_vertices);
    angular = sqrt(angular_error/num_vertices);
    ATE = sqrt(ATE/num_vertices);
end
