function plot_bigworld_with_distributedinfo_inverse_edge(graph, multi_graph, num_robots)

% Big world
figure('Position', [100, 100, 1200, 500]);
hold on;
x_big = [];
y_big = [];
for i = 1:length(graph.vertices_true)
    pose = graph.vertices_true{i};
    x_big(end+1) = pose(1);
    y_big(end+1) = pose(2);
end
plot(x_big, y_big, 'ko-', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', 'Big Map True Trajectory');

for edge_idx = 1:size(graph.edges, 1)
    v_i = graph.edges(edge_idx, 1);
    v_j = graph.edges(edge_idx, 2);
    if abs(v_j - v_i) > 1 
        x1 = graph.vertices_true{v_i}(1:2);
        x2 = graph.vertices_true{v_j}(1:2);
        plot([x1(1), x2(1)], [x1(2), x2(2)], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end
plot(x_big(1), y_big(1), 'go', 'MarkerSize', 10, 'LineWidth', 3, 'DisplayName', 'Start');

x_odom_big = [];
y_odom_big = [];
for i = 1:length(graph.vertices)
    pose_odom = graph.vertices{i};
    x_odom_big(end+1) = pose_odom(1);
    y_odom_big(end+1) = pose_odom(2);
end
plot(x_odom_big, y_odom_big, 'ro-', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', 'Big Map Odom Trajectory');

for edge_idx = 1:size(graph.edges, 1)
    v_i = graph.edges(edge_idx, 1);
    v_j = graph.edges(edge_idx, 2);
    if abs(v_j - v_i) > 1 
        x1_o = graph.vertices{v_i}(1:2);
        x2_o = graph.vertices{v_j}(1:2);
        plot([x1_o(1), x2_o(1)], [x1_o(2), x2_o(2)], 'r:', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

title('Original Square World (Multi-Robot)');
xlabel('X coordinate');
ylabel('Y coordinate');
legend;
grid on;
axis equal;
hold off;

% Multi robot world
figure('Position', [100, 100, 1200, 500]);
colors = ['r', 'g', 'b', 'm', 'c', 'k'];
hold on;

for r = 1 : num_robots
    x = [];
    y = [];
    for i = 1 : length(multi_graph(r).vertices_true)
        pose = multi_graph(r).vertices_true{i};
        x(end+1) = pose(1);
        y(end+1) = pose(2);
    end
    plot(x, y, [colors(r) 'o'], 'LineWidth',1.5,'MarkerSize',6, 'DisplayName', sprintf('Robot %d vertices', r));

    if ~isempty(x)
        plot(x(1), y(1), [colors(r) 'o'], 'MarkerSize', 20, 'LineWidth', 3, 'HandleVisibility', 'off');
    end

    intra = multi_graph(r).intra_edges;
    for intra_e = 1 : size(intra, 2)
        e_i = intra(1, intra_e);
        e_j = intra(2, intra_e);
        p_i = multi_graph(r).vertices_true{e_i}(1:2);
        p_j = multi_graph(r).vertices_true{e_j}(1:2);
        plot([p_i(1), p_j(1)], [p_i(2),p_j(2)], [colors(r) '-'], 'LineWidth', 1.2, 'HandleVisibility','off');
    end

    inter = multi_graph(r).inter_edges;
    for inter_e = 1 : size(inter, 2)
        inverse = inter(3, inter_e);
        if inverse == false
            e_i = inter(1, inter_e);
            lm_idx = inter(2, inter_e);

        else
            lm_idx = inter(1, inter_e);
            e_i = inter(2, inter_e);
        end
        p_i = multi_graph(r).vertices_true{e_i}(1:2);
        lm = multi_graph(r).lm_vertices_true{lm_idx}(1:2);
        plot([p_i(1), lm(1)], [p_i(2),lm(2)], [colors(r) ':'], 'LineWidth', 1.2, 'HandleVisibility','off');
    end

    lm_x = [];
    lm_y = [];
    for i = 1 : length(multi_graph(r).lm_vertices)
        lm_pose = multi_graph(r).lm_vertices_true{i};
        lm_x(end+1) = lm_pose(1);
        lm_y(end+1) = lm_pose(2);
    end
    plot(lm_x, lm_y, [colors(r) 's'], 'LineWidth',1.5,'MarkerSize',12, 'DisplayName', sprintf('Robot %d landmark', r));
end

title('Distributed Multi-Robot View');
xlabel('X coordinate');
ylabel('Y coordinate');
legend('Location', 'best');
grid on;
axis equal;
hold off


figure('Name','Intra Loop Closures per Robot','Position',[100 100 1200 300]);
ncols = min(3, num_robots);
nrows = ceil(num_robots / ncols);
tiledlayout(nrows, ncols, 'TileSpacing','compact','Padding','compact');

for r = 1:num_robots
    nexttile; hold on; grid on; axis equal;
    title(sprintf('Robot %d', r));
    xlabel('X'); ylabel('Y');

    xs = cellfun(@(p) p(1), multi_graph(r).vertices_true);
    ys = cellfun(@(p) p(2), multi_graph(r).vertices_true);
    plot(xs, ys, 'k.-', 'LineWidth', 1.5, 'DisplayName', 'trajectory');
    plot(xs(1), ys(1), 'go', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'start');
    plot(xs(end), ys(end), 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'end');

    E = multi_graph(r).intra_edges;
    lc_mask = abs(E(1,:) - E(2,:)) > 1;
    E_lc = E(:, lc_mask);
    for e = 1:size(E_lc, 2)
        i = E_lc(1, e); j = E_lc(2, e);
        plot([xs(i) xs(j)], [ys(i) ys(j)], 'r-', 'LineWidth', 1.5);
    end
end

% Show inter loop closures per robot 
figure('Name','Inter Loop Closures per Robot','Position',[100 100 1200 300]);
ncols = min(3, num_robots);
nrows = ceil(num_robots / ncols);
tiledlayout(nrows, ncols, 'TileSpacing','compact','Padding','compact');

for r = 1:num_robots
    nexttile; hold on; grid on; axis equal;
    title(sprintf('Robot %d (inter LCs)', r));
    xlabel('X'); ylabel('Y');

    xs = cellfun(@(p) p(1), multi_graph(r).vertices_true);
    ys = cellfun(@(p) p(2), multi_graph(r).vertices_true);
    plot(xs, ys, 'k.-', 'LineWidth', 1.5, 'DisplayName', 'trajectory');
    plot(xs(1), ys(1), 'go', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'start');
    plot(xs(end), ys(end), 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'end');
    
    Eint = multi_graph(r).inter_edges;              
    for e = 1:size(Eint, 2)
        inverse = Eint(3, e);
        if inverse == false
            i_loc  = Eint(1, e);
            lm_idx = Eint(2, e);
        else
            i_loc  = Eint(2, e);
            lm_idx = Eint(1, e);
        end
        
        p_local  = multi_graph(r).vertices_true{i_loc}(1:2);
        p_foreign = multi_graph(r).lm_vertices_true{lm_idx}(1:2);
        plot([p_local(1) p_foreign(1)], [p_local(2) p_foreign(2)], 'b--', 'LineWidth', 1.8, 'HandleVisibility','off');
    end

    if ~isempty(multi_graph(r).lm_vertices_true)
        lm_x = cellfun(@(p) p(1), multi_graph(r).lm_vertices_true);
        lm_y = cellfun(@(p) p(2), multi_graph(r).lm_vertices_true);
        plot(lm_x, lm_y, 'bs', 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'foreign landmarks');
    end
end




























