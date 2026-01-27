function [local_graph, F_now, norm_grad_F] = Optimization_RGN_linesearch(local_graph, alpha_init, grad_tol, lambda_reg)
% This is the main function for optimization action
    
    tau  = 1e-4;   
    beta = 0.5;
    found = false;
    l = 0;  
    alpha = alpha_init;  

    num_vertices = length(local_graph.vertices);
    anchor_first_vertex = isfield(local_graph,'anchor_first') && local_graph.anchor_first && isfield(local_graph,'robot_id') && local_graph.robot_id == 1;

    if anchor_first_vertex
        movable_indices = 5 : 4 * num_vertices;
    else
        movable_indices = 1 : 4 * num_vertices;
    end 

    [rgrad, rgnhess, F, ~] = rgn_gradhess_multi_J(local_graph);
    grad_norm = norm(rgrad(movable_indices));

    if grad_norm < grad_tol
        F_now = F;
        norm_grad_F = grad_norm;
        fprintf('[LineSearch] grad norm below threshold -> no update.\n');
        return;
    end

    H_k = rgnhess(movable_indices, movable_indices);
    g_k = rgrad(movable_indices);

    H_k_reg = H_k + lambda_reg * eye(size(H_k));
    s_k_trunc = - H_k_reg \ g_k;

    if anchor_first_vertex
        s_k = [zeros(4,1); s_k_trunc]; 
    else
        s_k = s_k_trunc;  
    end

    dd = dot(g_k, s_k_trunc);

    while ~found
        Exp_s_can = Exp_X_N(G_get_X(local_graph),beta^l * alpha * s_k);
        local_graph_candidate = G_set_X(local_graph, Exp_s_can);
        [~, ~, F_can, ~] = rgn_gradhess_multi_J(local_graph_candidate);

        if F_can <= F + tau * beta^l * alpha * dd
            alpha = beta^l * alpha;
            break;
        else
            l = l + 1;
            if beta^l * alpha < 1e-15
                warning('Armijo line search alpha is too small, breaking early.');
                break;
            end
        end
    end

    s_k_final = alpha * s_k;
    Exp_s = Exp_X_N(G_get_X(local_graph), s_k_final);
    local_graph = G_set_X(local_graph, Exp_s);

    [rgrad_final, ~, F_now, ~] = rgn_gradhess_multi_J(local_graph);
    norm_grad_F = norm(rgrad_final(movable_indices));
end








