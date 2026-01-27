function F = F_G_SE2(graph)
    F = 0.0;
    num_edges = length(graph.delta_poses_pudq);
    
    for ij=1:num_edges
        i = graph.edges(ij,1);
        j = graph.edges(ij,2);
        
        %Compute cost wrt the HTM SE2 metric
        T_ij = pose_to_SE2(graph.delta_poses{ij});
        T_i = pose_to_SE2(graph.vertices{i});
        T_j = pose_to_SE2(graph.vertices{j});

        eij = se2_vee(SE2_Log_1(SE2_inv(T_ij)*SE2_inv(T_i)*T_j));
        fij = eij'*graph.information_se2{ij}*eij;

        F = F + fij;
    end
    
    F = 0.5*F;
end
