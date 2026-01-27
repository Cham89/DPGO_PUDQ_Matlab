function G = G_set_X(G_0, X)
    num_vertices = length(G_0.vertices_pudq);
    
    G = G_0;
    for i=1:num_vertices
        i_index = 4*(i-1)+1;
        G.vertices_pudq{i} = X(i_index:i_index+3, 1);
        G.vertices{i} = pudq_to_pose(G.vertices_pudq{i});
    end
end

