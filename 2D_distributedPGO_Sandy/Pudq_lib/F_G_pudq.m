function F = F_G_pudq(graph)
    F = 0.0;
    num_edges = length(graph.delta_poses_pudq);
    
    %Loop through all edges and sum up their f_ij costs
    for ij=1:num_edges
        z_ij = graph.delta_poses_pudq{ij};
        x_i = graph.vertices_pudq{graph.edges(ij,1)};
        x_j = graph.vertices_pudq{graph.edges(ij,2)};

        eij = e_ij(z_ij, x_i, x_j);
        fij = eij' * graph.information_pudq{ij} * eij;
        F = F + fij;
    end
    
    F = 0.5*F;
end

