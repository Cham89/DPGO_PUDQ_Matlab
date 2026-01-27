function G = reorient_graph(G)
    num_vertices = length(G.vertices);
    
    x_0 = G.vertices_pudq{1};
    for i=1:num_vertices
        G.vertices_pudq{i} = pudq_mul(pudq_inv(x_0), G.vertices_pudq{i});
        
        % Update vertices to match
        G.vertices{i} = pudq_to_pose(G.vertices_pudq{i});
    end
end
