function [local_graph, F_now, norm_grad_F, Delta_out] = Optimization_RTR_inverse_edge(local_graph, Delta_in, grad_tol)

    rho_prime = 0.1;
    kappa = 0.1;
    theta = 1.0;
    max_inner = 4 * length(local_graph.vertices);
    Delta_bar = 1e6; 
    
    num_vertices = length(local_graph.vertices);
    anchor_first_vertex = isfield(local_graph,'anchor_first') && local_graph.anchor_first && isfield(local_graph,'robot_id') && local_graph.robot_id == 1;

    if anchor_first_vertex
        movable_indices = 5 : 4 * num_vertices;
    else
        movable_indices = 1 : 4 * num_vertices;
    end
    
    X_k = G_get_X(local_graph);
    [rgrad_k, gnhess_k, F_initial, P_X] = rgn_gradhess_multi_J_forinverse(local_graph);
    
    grad_norm = norm(rgrad_k(movable_indices));

    if grad_norm < grad_tol
        F_now = F_initial;
        norm_grad_F = grad_norm;
        Delta_out = Delta_in;
        fprintf(' Grad norm below threshold -> no update.\n');
        return;
    end

    H_k = gnhess_k(movable_indices, movable_indices);
    g_k = rgrad_k(movable_indices);
    P_X_trunc = P_X(movable_indices, movable_indices);

    [s_k_trunc, ~, ~, stop_tCG] = manopt_tCG(X_k, g_k, H_k, Delta_in, kappa, theta, P_X_trunc, max_inner);

    if anchor_first_vertex
        s_k = [zeros(4,1); s_k_trunc];
    else
        s_k = s_k_trunc;
    end

    pred_reduction = -(rgrad_k'*s_k + 0.5*s_k'*gnhess_k*s_k);
    X_candidate = Exp_X_N(X_k, s_k);
    graph_candidate = G_set_X(local_graph, X_candidate);
    [~, ~, F_candidate, ~] = rgn_gradhess_multi_J_forinverse(graph_candidate);

    actual_reduction = F_initial - F_candidate;

    if abs(pred_reduction) > 1e-9
        rho = actual_reduction / pred_reduction;
    else
        rho = (actual_reduction >= 0) * 1.0; 
    end

    if rho >= rho_prime
        local_graph = graph_candidate;
        F_now = F_candidate;
        [rgrad_final, ~, ~, ~] = rgn_gradhess_multi_J_forinverse(local_graph);
        norm_grad_F = norm(rgrad_final(movable_indices));
    else
        F_now = F_initial;
        norm_grad_F = grad_norm;
    end

    if (rho < 1/4)
        Delta_out = Delta_in / 4;
    elseif (rho > 3/4 && stop_tCG == 3) 
        Delta_out = min(2 * Delta_in, Delta_bar);
    else
        Delta_out = Delta_in;
    end
end