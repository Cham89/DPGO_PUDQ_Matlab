function Omega_ij = Omega_sparse_Multi(multi_graph)

    M = length(multi_graph.intra_info_pudq) + length(multi_graph.inter_info_pudq);
    num_triplets = 9*M;
    t_id = 0;
    i_triplets = zeros(num_triplets,1);
    j_triplets = zeros(num_triplets,1);
    v_triplets = zeros(num_triplets,1);
    all_info_pudq = [multi_graph.intra_info_pudq, multi_graph.inter_info_pudq];

    for ij=1:M
        info_pudq = all_info_pudq{ij};
        for i=1:3
            for j=1:3
                t_id = t_id+1;
                i_triplets(t_id) = 3*(ij-1)+i;
                j_triplets(t_id) = 3*(ij-1)+j;
                v_triplets(t_id) = info_pudq(i,j);
            end
        end
    end
    Omega_ij = sparse(i_triplets, j_triplets, v_triplets, 3*M, 3*M);
end
