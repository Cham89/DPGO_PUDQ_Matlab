function RPE = RPE_G_pudq(graph)
    RPE = 0.0;
    sz_edges = size(graph.delta_poses_pudq);
    num_edges = sz_edges(2);
    for i=1:num_edges
        z_ij_true = graph.delta_poses_pudq_true{i};

        x_i = graph.vertices_pudq{graph.edges(i,1)};
        x_j = graph.vertices_pudq{graph.edges(i,2)};

        eij = e_ij(z_ij_true, x_i, x_j);

        RPE = RPE + eij' * eij;
    end
        RPE = sqrt(RPE/num_edges);
end
