%Absolute trajectory error
function ATE = ATE_G_pudq(graph)
    ATE = 0.0;
    num_vertices = length(graph.vertices);
    for i=1:num_vertices
        x_i_true = graph.vertices_pudq_true{i};
        x_i_hat = graph.vertices_pudq{i};
        id = pudq_identity();

        eij = e_ij(x_i_true, id, x_i_hat);

        ATE = ATE + eij' * eij;
    end
       
    ATE = sqrt(ATE/num_vertices);
end
