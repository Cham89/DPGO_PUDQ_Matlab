function boumal_check_multirobot(local_graph_r, robot_id)
    fprintf('\n=== BOUMAL CHECK FOR ROBOT %d ===\n', robot_id);
    num_vertices = length(local_graph_r.vertices);
    anchor_first_vertex = isfield(local_graph_r,'anchor_first') && local_graph_r.anchor_first && isfield(local_graph_r,'robot_id') && local_graph_r.robot_id == 1;
    
    if anchor_first_vertex
        movable_indices = 5 : 4 * num_vertices;
    else
        movable_indices = 1 : 4 * num_vertices;
    end

    X_current = G_get_X(local_graph_r);
    [rgradF, ~, F_current, P_X] = rgn_gradhess_multi_J_forinverse(local_graph_r);

    rgrad_movable = rgradF(movable_indices);
    
    fprintf('Current cost: %.6e\n', F_current);
    fprintf('Gradient norm: %.6e\n', norm(rgrad_movable));
    
    V_movable = randn(length(movable_indices), 1);
    V_movable = V_movable / norm(V_movable);
    
    if anchor_first_vertex
        V_full = [zeros(4,1); V_movable];  
    else
        V_full = V_movable;
    end
    
    N_t = 500;
    t = logspace(-8, 0, N_t);
    E_g = zeros(1, N_t);
    fprintf('Testing %d step sizes...\n', N_t);
    
    valid_count = 0;
    for j = 1:N_t
        X_new = Exp_X_N(X_current, t(j) * V_full);
        local_graph_test = G_set_X(local_graph_r, X_new);
        [~, ~, F_new] = rgn_gradhess_multi_J_forinverse(local_graph_test);
        
        predicted_change = t(j) * dot(rgrad_movable, V_movable);
        actual_change = F_new - F_current;
        E_g(j) = abs(actual_change - predicted_change);
        valid_count = valid_count + 1;
    end
    
    valid_mask = isfinite(E_g) & E_g > 0;
    if sum(valid_mask) < 20
        fprintf('ERROR: Only %d valid points', sum(valid_mask));
        return;
    end
    
    valid_t = t(valid_mask);
    valid_E_g = E_g(valid_mask);
    
    mid_start = max(1, floor(length(valid_t)/3));
    mid_end = min(length(valid_t), floor(2*length(valid_t)/3));
    fit_t = valid_t(mid_start:mid_end);
    fit_E_g = valid_E_g(mid_start:mid_end);
    
    log_t = log10(fit_t);
    log_E = log10(fit_E_g);
    p = polyfit(log_t, log_E, 1);
    slope = p(1);
    
    slope_2 = polyval([2.0, 0.0], log10(t));
    log_slope_2 = 10.^(slope_2);
    
    figure
    hold on
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    
    plot_title = ['Riemannian Gradient Check: ' num2str(num_vertices) ' vertices'];
    title(plot_title)
    grid on
    box on
    axis equal
    
    plot(t, E_g, 'DisplayName', 'E_g', 'LineWidth', 1.5)
    plot(t, log_slope_2, '-r', 'Displayname', 'Slope 2')
    legend('Location', 'southeast')
    
    hold off
    
    fprintf('\n=== RESULTS ===\n');
    fprintf('Computed slope: %.3f\n', slope);
    fprintf('Valid points: %d/%d\n', sum(valid_mask), N_t);
    
    if abs(slope - 2.0) < 0.4
        fprintf('✓ GRADIENT IS CORRECT (slope ≈ 2.0)\n');
    else
        fprintf(' Blaaa ');
    end
    
end