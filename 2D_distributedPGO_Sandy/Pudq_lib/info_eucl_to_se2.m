function info_mat_se2 = info_eucl_to_se2(info_mat_eucl, theta)
    M = [sinc1(theta), -cosc(theta), 0; cosc(theta), sinc1(theta), 0; 0, 0, 1];
    info_mat_se2 = M'*info_mat_eucl*M;
    info_mat_se2 = tril(info_mat_se2)+tril(info_mat_se2,-1)';
end
