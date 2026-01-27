function M = create_block_jacobi_preconditioner(H)

    
    n_dim = size(H, 1);
    if mod(n_dim, 4) ~= 0
        error('Hessian 矩陣的維度必須是 4 的倍數。');
    end
    num_blocks = n_dim / 4;
    
    
    num_triplets = 16 * num_blocks;
    I = zeros(num_triplets, 1);
    J = zeros(num_triplets, 1);
    S = zeros(num_triplets, 1);
    
    
    [ii, jj] = meshgrid(1:4, 1:4);
    
    
    for i = 1:num_blocks
        
        idx = (4*(i-1)+1) : (4*i);
        
        
        block = full(H(idx, idx));
        
        
        inv_block = inv(block + eye(4) * 1e-9);
        
        
        triplet_indices = (16*(i-1)+1) : (16*i);
        
        
        I(triplet_indices) = idx(ii(:));
        J(triplet_indices) = idx(jj(:));
        
       
        S(triplet_indices) = inv_block(:);
    end
    
   
    M = sparse(I, J, S, n_dim, n_dim);
end