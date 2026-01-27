function [local_graph, F_now, norm_grad_F] = Optimization_RLM_inverse_edge(local_graph, mu_in, grad_tol)
    eta = 0.2;
    beta = 5;
    max_lm_trials = 300;

    num_vertices = length(local_graph.vertices);
    anchor_first_vertex = isfield(local_graph,'anchor_first') && local_graph.anchor_first && isfield(local_graph,'robot_id') && local_graph.robot_id == 1;
    if anchor_first_vertex
        movable_indices = 5 : 4 * num_vertices;
    else
        movable_indices = 1 : 4 * num_vertices;
    end

    [rgrad, rgnhess, F_initial, ~] = rgn_gradhess_multi_J_forinverse(local_graph);
    grad_norm = norm(rgrad(movable_indices));

    if grad_norm < grad_tol
        F_now = F_initial;
        norm_grad_F = grad_norm;  
        fprintf(' Grad norm below threshold -> no update.\n');
        return;
    end

    H_k = rgnhess(movable_indices, movable_indices);
    g_k = rgrad(movable_indices);
    normF2 = 2.0 * F_initial;
     mu_trial = mu_in;
    %lambda = mu_in;
    lambda = mu_trial * normF2;

    for lm_trial = 1:max_lm_trials
        
      
        H_k_reg = H_k + lambda * eye(size(H_k));
        % H_k_reg_full = rgnhess + lambda * eye(size(rgnhess));
        s_k_trunc = - H_k_reg \ g_k;
        
        if anchor_first_vertex
            s_k = [zeros(4,1); s_k_trunc];
        else
            s_k = s_k_trunc;
        end
        
        X_candidate = Exp_X_N(G_get_X(local_graph), s_k);
        graph_candidate = G_set_X(local_graph, X_candidate);
        [~, ~, F_candidate, ~] = rgn_gradhess_multi_J_forinverse(graph_candidate);
        % theta_k_0 = F_initial^2;
        % theta_k_s = F_initial^2 + 2*rgrad'*s_k + s_k'*rgnhess*s_k + lambda*(s_k'*s_k);
        actual_reduction = F_initial - F_candidate;
        %pred_reduction =  0.5*(theta_k_0 - theta_k_s);
        pred_reduction = -dot(g_k, s_k_trunc) - 0.5 * dot(s_k_trunc, H_k_reg * s_k_trunc);
        
        if pred_reduction > 0
            rho = actual_reduction / pred_reduction;
        else
            rho = -inf;
        end
        if rho >= eta
            local_graph = graph_candidate;
            F_now = F_candidate;
            [rgrad_final, ~, ~, ~] = rgn_gradhess_multi_J_forinverse(local_graph);
            norm_grad_F = norm(rgrad_final(movable_indices));
            return;
        else
            mu_trial = mu_trial * beta;
            %lambda = lambda * beta;
        end
    end
    

    warning('RLM failed to find a successful step after %d trials. Lambda is now very large.', max_lm_trials);
    F_now = F_initial;
    norm_grad_F = grad_norm;
end












