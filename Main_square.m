% Main function for the square world test
clear; clc; close all

pudq_pgo_lib_path = [pwd '/Pudq_lib'];
addpath(pudq_pgo_lib_path)

% Parameters
seed             = 22316;
rng(seed);

inverse = true;
return_early = false;

num_robots       = 3;
sigma_pudq       = 1e-2;
side_len         = 5;
gap              = 2;
inter_lc_prob    = 0.05;
lambda            = 1;
N_max            = 800000;
alpha_init       = 1e-2;
grad_tol         = 1e-2;
lambda_reg       = 1e-4;

graph = generate_pudq_squareworld_multirobot_2d(num_robots, sigma_pudq, side_len , gap, inter_lc_prob);

disp('Square world graph generation complete.');

if inverse 
    multi_graph = distributebiggraph_withinverseedge(graph, num_robots);
else
    multi_graph = distributebiggraph(graph, num_robots);
end

if inverse 
    plot_bigworld_with_distributedinfo_inverse_edge (graph, multi_graph, num_robots);
    if return_early
        return;
    end
else
    plot_bigworld_with_distributedinfo(graph, multi_graph, num_robots);
    if return_early
        return;
    end
end

robot = struct;
for r = 1: num_robots
    local_graph = build_local_subgraph(multi_graph, r);

    if inverse 
        [rgrad_0, ~, F_0] = rgn_gradhess_multi_J_forinverse(local_graph);
    else
        [rgrad_0, ~, F_0] = rgn_gradhess_multi_J(local_graph);
    end
    gradnorm_0 = norm(rgrad_0);

    robot(r).id = r;
    robot(r).local_graph = local_graph;
    robot(r).rgn_cost = F_0;
    robot(r).rgn_gradnorm = gradnorm_0;
    robot(r).converge = false;
end


global_time = 0;
robot_time = exprnd (1/lambda, [1, num_robots]);

for iter = 1 : N_max
    [global_time, r] = min(robot_time);

    if inverse 
        [local_graph, F_now, norm_grad_F] = Optimization_RGN_linesearch_inverse_edge(robot(r).local_graph, alpha_init, grad_tol, lambda_reg);
    else
        [local_graph, F_now, norm_grad_F] = Optimization_RGN_linesearch(robot(r).local_graph, alpha_init, grad_tol, lambda_reg);
    end
        
    robot(r).local_graph      = local_graph;
    robot(r).rgn_cost(end+1)     = F_now;
    robot(r).rgn_gradnorm(end+1) = norm_grad_F;

    robot = shareupdateinfo(r, local_graph, robot);

    robot_time(r) = global_time + exprnd(1/lambda);
    fprintf('Iter=%3d -> Robot %d updated  (F=%.3g, grad=%.3g)\n', iter, r, F_now, norm_grad_F);

    
    if inverse 
        robot = check_all_robotgrad_inverse_edge(robot, grad_tol);
    else
        robot = check_all_robotgrad(robot, grad_tol);
    end

    if all([robot.converge])
        fprintf('All %d robots converged -> terminating at iter=%d.\n', numel(robot), iter);
        break;
    end

end

fprintf('\n*** Done with %d asynchronous update events! ***\n', N_max);


figure('Name','Per-Robot Cost','NumberTitle','off');
hold on; grid on;
colors = lines(num_robots);

for r = 1:num_robots
    N = length(robot(r).rgn_cost) - 1;
    plot_indices = linspace(0, N, N+1);
    
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
    plot_indices = linspace(0, N, N+1);
    
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
% 
% figure('Name','Per-Robot Grad Norm (Global)','NumberTitle','off');
% set(gca, 'YScale', 'log');

% hold on; grid on;
% for r = 1:num_robots
%     semilogy(1:numel(robot(r).rgn_gradnorm_global), robot(r).rgn_gradnorm_global, '-', 'LineWidth', 1.5, ...
%              'Color', colors(r,:), 'DisplayName', sprintf('Robot %d', r));  
% end
% xlabel('Global async event'); ylabel('RGN gradient norm');
% legend('show', 'Location','best');
% ax = gca;
% title('Per-robot gradient (every iteration)');
% hold off;














