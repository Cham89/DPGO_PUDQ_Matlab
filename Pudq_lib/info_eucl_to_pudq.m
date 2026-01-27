function info_mat_pudq = info_eucl_to_pudq(info_mat_eucl, theta_error)

    alpha = theta_error/2.0;
    beta = cos(alpha)/sinc1(alpha);
    
    M = [1, 0, 0; 0, beta, alpha; 0, -alpha, beta]; %PUDQ transformation matrix
    B = [0, 0, 1; 1, 0, 0; 0, 1, 0];
    info_mat_pudq = 4*B*((M' \ info_mat_eucl) / M)*B';
    info_mat_pudq = tril(info_mat_pudq)+tril(info_mat_pudq,-1)';
end
