function info_mat_se2 = info_pudq_to_se2(info_mat_pudq)
    
    %Coordinate remapping matrix
    B = [0, 0, 1; 1, 0, 0; 0, 1, 0];
    
    %Perform coordinate remapping, then transform to pudq space
    info_mat_se2 = 1/4*B'*info_mat_pudq*B;

    %Use the lower triangular part to ensure symmetricity
%     info_mat_eucl = tril(info_mat_eucl)+tril(info_mat_eucl,-1)';
end
