function X = G_get_X(G)
    num_vertices = length(G.vertices);
    X = zeros(num_vertices*4, 1);

    for i=1:num_vertices
        i_index = 4*(i-1)+1;
        X(i_index:i_index+3,1) = G.vertices_pudq{i};
    end
end
