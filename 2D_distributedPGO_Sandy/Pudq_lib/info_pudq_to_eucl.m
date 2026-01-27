function info_mat_eucl = info_pudq_to_eucl(info_mat_pudq, theta_error)

    alpha = theta_error/2.0;
    beta = cos(alpha)/sinc1(alpha);
    
    %PUDQ transformation matrix
    M = [1, 0, 0; 0, beta, alpha; 0, -alpha, beta];
    
    %Coordinate remapping matrix
    B = [0, 0, 1; 1, 0, 0; 0, 1, 0];
    
    %Perform coordinate remapping, then transform to eucl space
    info_mat_eucl = 0.25*M'*B'*info_mat_pudq*B*M;

    %Use the lower triangular part to ensure symmetricity
    % info_mat_eucl = tril(info_mat_eucl)+tril(info_mat_eucl,-1)';
end