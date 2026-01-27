clear, clc, close all;

run('Cartan-Sync\MATLAB\setup.m');

pudq_pgo_lib_path = [pwd '/Pudq_lib'];
addpath(pudq_pgo_lib_path);

riem_opt_path = [pwd '/riem_opt_lib'];
addpath(riem_opt_path);

se2_lib_path = [pwd '/se2_lib'];
addpath(se2_lib_path);

Manopt_opts = {};
SE_Sync_opts.init = 'chordal';
SE_Sync_opts.Precon = true;

inverse = true;
return_early = false;  
num_robots = 3;
lambda                 = 1;
N_max                  = 800000;
alpha_init             = 1e-1;
grad_tol               = 1e-2;
lambda_reg             = 1e-2;


gt_g2o_file = 'C:\Users\tcche\OneDrive\桌面\2D_distributedPGO_Sandy\2D_distributedPGO_Sandy\data\Grid1000_ground_truth.g2o';
GT_graph = load_g2o_file_pudq_gt (gt_g2o_file);

odom_g2o_file = 'C:\Users\tcche\OneDrive\桌面\2D_distributedPGO_Sandy\2D_distributedPGO_Sandy\data\Grid1000_1.g2o';
odom_graph = load_g2o_file_pudq (odom_g2o_file);

graph = combine_gt_odom_g2o(odom_graph, GT_graph);

graph_RGN = graph;
num_vertices = length(graph.vertices);
[Rchordal, tchordal] = chordal_initialization(graph_RGN);
for i=1:num_vertices
    R_i_index = 2*(i-1)+1;
    R_i = Rchordal(:,R_i_index:R_i_index+1);
    t_i = tchordal(:,i);

    theta_i = theta_to_R(R_i);
    pose_i = [t_i; theta_i];
    pudq_i = pose_to_pudq(pose_i);

    G_rtr.vertices{i} = pose_i;
    G_rtr.vertices_pudq{i} = pudq_i;
end

if inverse 
    multi_graph_RGN = distributebiggraph_withinverseedge(graph_RGN, num_robots);
else
    multi_graph_RGN = distributebiggraph(graph_RGN, num_robots);
end

if inverse 
    plot_bigworld_with_distributedinfo_inverse_edge(graph_RGN, multi_graph_RGN, num_robots);
    if return_early
        return;
    end
else
    plot_bigworld_with_distributedinfo(graph_RGN, multi_graph_RGN, num_robots);
    if return_early
        return;
    end
end


robot = struct;
for r = 1: num_robots
    local_graph = build_local_subgraph(multi_graph_RGN, r);

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
iter = 0;


while true
    
    iter = iter + 1;
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

