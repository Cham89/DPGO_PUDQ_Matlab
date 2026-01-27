function plot_graph_trajectory_with_ground_truth(graph, title_str, show_ground_truth)
    if nargin < 3
        show_ground_truth = false;
    end
    
    hold on;
    grid on;
    axis equal;
    
    % Extract coordinates from estimated vertices
    num_vertices = length(graph.vertices);
    x_coords = zeros(1, num_vertices);
    y_coords = zeros(1, num_vertices);
    
    for i = 1:num_vertices
        pose = graph.vertices{i};
        x_coords(i) = pose(1);
        y_coords(i) = pose(2);
    end
    
    % Show ground truth trajectory if available and requested
    if show_ground_truth && isfield(graph, 'vertices_true')
        x_true = zeros(1, num_vertices);
        y_true = zeros(1, num_vertices);
        
        for i = 1:num_vertices
            pose_true = graph.vertices_true{i};
            x_true(i) = pose_true(1);
            y_true(i) = pose_true(2);
        end
        
        % Plot ground truth
        plot(x_true, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'Ground Truth');
        
        % Plot ground truth start/end
        plot(x_true(1), y_true(1), 'ko', 'MarkerSize', 8, 'LineWidth', 2, ...
             'DisplayName', 'True Start', 'MarkerFaceColor', 'k');
        plot(x_true(end), y_true(end), 'ks', 'MarkerSize', 8, 'LineWidth', 2, ...
             'DisplayName', 'True End', 'MarkerFaceColor', 'k');
    end
    
    % Plot the estimated trajectory
    plot(x_coords, y_coords, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3, ...
         'DisplayName', 'Estimated', 'MarkerFaceColor', 'b');
    
    % Plot estimated start/end points
    plot(x_coords(1), y_coords(1), 'go', 'MarkerSize', 8, 'LineWidth', 2, ...
         'DisplayName', 'Est. Start', 'MarkerFaceColor', 'g');
    plot(x_coords(end), y_coords(end), 'ro', 'MarkerSize', 8, 'LineWidth', 2, ...
         'DisplayName', 'Est. End', 'MarkerFaceColor', 'r');
    
    % Plot loop closure edges
    if isfield(graph, 'edges')
        for edge_idx = 1:size(graph.edges, 1)
            v1 = graph.edges(edge_idx, 1);
            v2 = graph.edges(edge_idx, 2);
            
            % Only show non-consecutive edges (loop closures)
            if abs(v1 - v2) > 1
                p1 = [x_coords(v1), y_coords(v1)];
                p2 = [x_coords(v2), y_coords(v2)];
                plot([p1(1), p2(1)], [p1(2), p2(2)], 'r--', 'LineWidth', 1, ...
                     'HandleVisibility', 'off', 'Color', [0.8, 0.2, 0.2, 0.6]);
            end
        end
    end
    
    % Add title and labels
    title(title_str, 'FontSize', 12);
    xlabel('X coordinate');
    ylabel('Y coordinate');
    
    % Add legend
    legend('show', 'Location', 'best');
    
    hold off;
end
