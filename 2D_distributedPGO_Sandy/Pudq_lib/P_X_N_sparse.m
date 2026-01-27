function P = P_X_N_sparse(X)

    N = length(X);
    num_triplets = 16*N;
    t_id = 0;
    i_triplets = zeros(num_triplets, 1);
    j_triplets = zeros(num_triplets, 1);
    v_triplets = zeros(num_triplets, 1);


    for v=1:N
        x_i = X{v};
        P_x = [eye(2) - x_i(1:2)*x_i(1:2)', zeros(2,2); zeros(2,2), eye(2)];

        for i=1:4
            for j=1:4
                t_id = t_id+1;
                i_triplets(t_id) = 4*(v-1)+i;
                j_triplets(t_id) = 4*(v-1)+j;
                v_triplets(t_id) = P_x(i,j);
            end
        end

    end
    
    P = sparse(i_triplets, j_triplets, v_triplets, 4*N, 4*N);
end
