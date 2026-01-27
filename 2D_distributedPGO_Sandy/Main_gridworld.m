% Main function for the square world test
clear; clc; close all

pudq_pgo_lib_path = [pwd '/Pudq_lib'];
addpath(pudq_pgo_lib_path);

% Parameters
seed             = 224832901;
rng(seed);

iter = 0;
inverse = true;
return_early = false;  
num_robots             = 3;
num_vertices_per_robot = 10;
lc_dist                = 3.0;
lc_prob                = 0.3;
sigma_pudq             = 1e-3;
lambda                 = 1;
N_max                  = 800000;
alpha_init             = 1e-2;
grad_tol               = 1e-2;
lambda_reg             = 1e-4;

graph = generate_pudq_gridworld_multirobot_2d(num_vertices_per_robot, num_robots, sigma_pudq, lc_dist, lc_prob);

disp('Square world graph generation complete.');

multi_graph = distributebiggraph_withinverseedge(graph, num_robots);

plot_bigworld_with_distributedinfo_inverse_edge(graph, multi_graph, num_robots);
if return_early
    return;
end



robot = struct;
for r = 1: num_robots
    local_graph = build_local_subgraph(multi_graph, r);
    [rgrad_0, ~, F_0] = rgn_gradhess_multi_J_forinverse(local_graph);

    gradnorm_0 = norm(rgrad_0);

    robot(r).id = r;
    robot(r).local_graph = local_graph;
    robot(r).rgn_cost = F_0;
    robot(r).rgn_gradnorm = gradnorm_0;
    robot(r).converge = false;
end

global_time = 0;
robot_time = exprnd (1/lambda, [1, num_robots]);

while true
    iter = iter + 1;
    
    [global_time, r] = min(robot_time);

    [local_graph, F_now, norm_grad_F] = Optimization_RGN_linesearch_inverse_edge(robot(r).local_graph, alpha_init, grad_tol, lambda_reg);

    robot(r).local_graph      = local_graph;
    robot(r).rgn_cost(end+1)     = F_now;
    robot(r).rgn_gradnorm(end+1) = norm_grad_F;

    robot = shareupdateinfo(r, local_graph, robot);

    robot_time(r) = global_time + exprnd(1/lambda);
    fprintf('Iter=%3d -> Robot %d updated  (F=%.3g, grad=%.3g)\n', iter, r, F_now, norm_grad_F);
 
    robot = check_all_robotgrad_inverse_edge(robot, grad_tol);
    if all([robot.converge])
        fprintf('All %d robots converged -> terminating at iter=%d.\n', numel(robot), iter);
        break;
    end

end

fprintf('\n*** Done with %d asynchronous update events! ***\n', N_max);
graph_optimized = reconstructglobalgraph(robot, graph);
graph_multirobot = reorient_graph(graph_optimized);
figure('Name', 'Detailed Comparison', 'Position', [100, 100, 1800, 600]);
plot_graph_trajectory_with_ground_truth(graph_multirobot, 'Distributed Multi-Robot', true);
sgtitle('Results vs Ground Truth', 'FontSize', 16);


figure('Name','Per-Robot Cost','NumberTitle','off');
hold on; grid on;
colors = lines(num_robots);

for r = 1:num_robots
    N = length(robot(r).rgn_cost) - 1;
    plot_indices = 1:length(robot(r).rgn_cost);
    
    plot( plot_indices, robot(r).rgn_cost,'Color', colors(r,:), 'LineWidth', 1.25, 'DisplayName', sprintf('Robot %d', r) );
end

% set(gca, 'YScale', 'log');
title('Per-Robot Cost vs. iteration (Asynchronous Poisson updates)', 'FontSize', 14);
xlabel('Local iteration count', 'FontSize', 12);
ylabel('RGD Cost',              'FontSize', 12);
legend('show', 'Location','best');
ax = gca;
ax.YAxis.TickLabelFormat = '10^{%g}';

hold off;

figure('Name','Per-Robot Grad Norm','NumberTitle','off');
hold on; grid on;

for r = 1:num_robots
    N = length(robot(r).rgn_gradnorm) - 1;
    plot_indices = 1:length(robot(r).rgn_gradnorm);
    
    plot( plot_indices, robot(r).rgn_gradnorm, 'Color', colors(r,:),'LineWidth', 1.25, 'DisplayName', sprintf('Robot %d', r) );
end

set(gca, 'YScale', 'log');
title('Per-Robot Gradient Norm (Asynchronous Poisson)', 'FontSize', 14);
xlabel('Local iteration count', 'FontSize', 12);
ylabel('RGD Gradient norm',    'FontSize', 12);
legend('show', 'Location','best');
ax = gca;
ax.YAxis.TickLabelFormat = '10^{%g}';

hold off;
