clear, clc, close all;

run('/home/hsuanpin/Desktop/2D_distributedPGO_Sandy/Cartan-Sync/MATLAB/setup.m');

pudq_pgo_lib_path      = [pwd '/Pudq_lib'];
addpath(pudq_pgo_lib_path);

riem_opt_path          = [pwd '/riem_opt_lib'];
addpath(riem_opt_path);

se2_lib_path           = [pwd '/se2_lib'];
addpath(se2_lib_path);

Manopt_opts            = {};
SE_Sync_opts.init      = 'chordal';
SE_Sync_opts.Precon    = true;


% Parameters
seed                   = 224832901;
rng(seed);

iter                   = 0;
return_early           = false;  
num_robots             = 2;
num_vertices_per_robot = 1750;
lc_dist                = 3.0;
lc_prob                = 0.02;
sigma_pudq             = 1e-3;
lambda                 = 1;
alpha_init             = 1e-1;
grad_tol               = 1e-2;
lambda_reg             = 1e-3;

graph = generate_pudq_gridworld_multirobot_2d(num_vertices_per_robot, num_robots, sigma_pudq, lc_dist, lc_prob);
disp('Square world graph generation complete.');

%%%%%%%%%% Chordal Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = graph;
num_vertices = length(graph.vertices);
[Rchordal, tchordal] = chordal_initialization(G);
for i=1:num_vertices
    R_i_index = 2*(i-1)+1;
    R_i = Rchordal(:,R_i_index:R_i_index+1);
    t_i = tchordal(:,i);

    theta_i = theta_to_R(R_i);
    pose_i = [t_i; theta_i];
    pudq_i = pose_to_pudq(pose_i);

    G.vertices{i} = pose_i;
    G.vertices_pudq{i} = pudq_i;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

multi_graph = distributebiggraph_withinverseedge(G, num_robots);
% plot_bigworld_with_distributedinfo_inverse_edge(G, multi_graph, num_robots);
if return_early
   return;
end

%%%%%%%%%% Cartan-Sync %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cartan_problem = SE_Sync_Cartan(G, Manopt_opts, SE_Sync_opts);
[SDPval_cartan, Yopt_cartan, Xhat_cartan, Fxhat_cartan, cartan_info] = run(cartan_problem);

%Convert SE sync data to graph
G_Cartan = se_sync_to_graph(G, Xhat_cartan);

%Get the final SE-Sync cost
F_Cartan = F_G_pudq(G_Cartan);

%Get SE-Sync time to convergence
cartan_time = sum([cartan_info.mat_construct_times, cartan_info.init_time, cartan_info.optimization_times, cartan_info.min_eig_vals, cartan_info.min_eig_times]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% DPGO RGN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    %[local_graph, F_now, norm_grad_F] = Optimization_RGN_linesearch_inverse_edge(robot(r).local_graph, alpha_init, grad_tol, lambda_reg);
    [local_graph, F_now, norm_grad_F] = Optimization_RLM_inverse_edge(robot(r).local_graph, lambda_reg, grad_tol);

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

graph_optimized = reconstructglobalgraph(robot, G);
fprintf('\n*** Done with asynchronous update events! ***\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G_multirobot = reorient_graph(graph_optimized);
G_Cartan         = reorient_graph(G_Cartan);

[original_rgrad, ~, original_F, ~] = rgn_gradhess_J(G);
original_se2_cost = F_G_SE2(G);
[multirobot_rgrad, ~, multirobot_F, ~] = rgn_gradhess_J(G_multirobot);
multirobot_se2_cost = F_G_SE2(G_multirobot);
[cartan_rgrad, ~, cartan_F, ~] = rgn_gradhess_J(G_Cartan);
cartan_se2_cost = F_G_SE2(G_Cartan);

multirobot_rpe_pudq = RPE_G_pudq(G_multirobot);
cartan_rpe_pudq = RPE_G_pudq(G_Cartan);

[multirobot_rpe_eucl, multirobot_rpe_linear, multirobot_rpe_angular] = RPE_G_eucl(G_multirobot);
[cartan_rpe_eucl, cartan_rpe_linear, cartan_rpe_angular] = RPE_G_eucl(G_Cartan);

multirobot_ate_pudq = ATE_G_pudq(G_multirobot);
cartan_ate_pudq = ATE_G_pudq(G_Cartan);

[multirobot_ate_eucl, multirobot_ate_linear, multirobot_ate_angular] = ATE_G_eucl(G_multirobot);
[cartan_ate_eucl, cartan_ate_linear, cartan_ate_angular] = ATE_G_eucl(G_Cartan);

fprintf('Multi-Robot Results:\n');
fprintf('  SE2 Cost: %.6f\n', multirobot_se2_cost);
fprintf('  PUDQ Cost: %.6f\n', multirobot_F);
fprintf('  Grad: %.6f\n', norm(multirobot_rgrad));
fprintf('  RPE PUDQ: %.6f\n', multirobot_rpe_pudq);
fprintf('  RPE Euclidean: %.6f (Linear: %.6f, Angular: %.6f)\n', multirobot_rpe_eucl, multirobot_rpe_linear, multirobot_rpe_angular);
fprintf('  ATE PUDQ: %.6f\n', multirobot_ate_pudq);
fprintf('  ATE Euclidean: %.6f (Linear: %.6f, Angular: %.6f)\n', multirobot_ate_eucl, multirobot_ate_linear, multirobot_ate_angular);

fprintf('Cartin-sync Results:\n');
fprintf('  SE2 Cost: %.6f\n', cartan_se2_cost);
fprintf('  PUDQ Cost: %.6f\n', cartan_F);
fprintf('  Grad: %.6f\n', norm(cartan_rgrad));
fprintf('  RPE PUDQ: %.6f\n', cartan_rpe_pudq);
fprintf('  RPE Euclidean: %.6f (Linear: %.6f, Angular: %.6f)\n', cartan_rpe_eucl, cartan_rpe_linear, cartan_rpe_angular);
fprintf('  ATE PUDQ: %.6f\n', cartan_ate_pudq);
fprintf('  ATE Euclidean: %.6f (Linear: %.6f, Angular: %.6f)\n', cartan_ate_eucl, cartan_ate_linear, cartan_ate_angular);


% Compare with GT
figure('Name', 'Detailed Comparison', 'Position', [100, 100, 1800, 600]);
subplot(1, 3, 1);
plot_graph_trajectory_with_ground_truth(G, 'Original (Odometry)', true);
subplot(1, 3, 2);
plot_graph_trajectory_with_ground_truth(G_Cartan, 'Cartan-Sync', true);
subplot(1, 3, 3);
plot_graph_trajectory_with_ground_truth(G_multirobot, 'Distributed Multi-Robot', true);
sgtitle('Results vs Ground Truth', 'FontSize', 16);



