% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT: run_analysis.m
% Purpose: Load a converged simulation result and run
%          the final plotting and analysis against
%          the original graph and Cartan-Sync.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all;

fprintf('*** RUNNING FINAL ANALYSIS ***\n');

% --- 1. SETUP ---
% Add the paths to your libraries
run('Cartan-Sync\MATLAB\setup.m');
addpath([pwd '/Pudq_lib']);
addpath([pwd '/riem_opt_lib']);
addpath([pwd '/se2_lib']);

% Define parameters needed for setup
num_robots = 2;
Manopt_opts = {};
SE_Sync_opts.init = 'chordal';
SE_Sync_opts.Precon = true;

% --- 2. LOAD ORIGINAL DATA ---
% We must re-load the original data to get 'G' for comparison
fprintf('Loading g2o files...\n');
gt_g2o_file = 'C:\Users\hchen804\OneDrive - Georgia Institute of Technology\Desktop\2D_distributedPGO_Sandy\data\Grid1000_ground_truth.g2o';
odom_g2o_file = 'C:\Users\hchen804\OneDrive - Georgia Institute of Technology\Desktop\2D_distributedPGO_Sandy\data\Grid1000_1.g2o';
GT_graph = load_g2o_file_pudq_gt (gt_g2o_file);
odom_graph = load_g2o_file_pudq (odom_g2o_file);
graph = combine_gt_odom_g2o(odom_graph, GT_graph);

% --- 3. RE-RUN CHORDAL INITIALIZATION (to get G) ---
fprintf('Running Chordal Initialization...\n');
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

% --- 4. RE-RUN CARTAN-SYNC (to get G_Cartan) ---
fprintf('Running Cartan-Sync for comparison...\n');
cartan_problem = SE_Sync_Cartan(G, Manopt_opts, SE_Sync_opts);
[~, ~, Xhat_cartan, ~, ~] = run(cartan_problem);
G_Cartan = se_sync_to_graph(G, Xhat_cartan);
F_Cartan = F_G_pudq(G_Cartan); % Used in analysis

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- 5. LOAD YOUR CONVERGED RESULTS ---
% This loads the 'robot' struct that your
% 'resume_simulation.m' script saved.
fprintf('Loading converged robot state from file...\n');
load('converged_robot_state.mat');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- 6. RUN FINAL PLOTTING AND ANALYSIS ---
% This is the *exact* code from the end of your main script.

fprintf('Generating plots and calculating final results...\n');

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
ylabel('RGD Cost', 'FontSize', 12);
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
ylabel('RGD Gradient norm', 'FontSize', 12);
legend('show', 'Location','best');
ax = gca;
ax.YAxis.TickLabelFormat = '10^{%g}';
hold off;

graph_optimized = reconstructglobalgraph(robot, G);
fprintf('\n*** Done with asynchronous update events! ***\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G_multirobot = reorient_graph(graph_optimized);
G_Cartan = reorient_graph(G_Cartan);

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
fprintf('  RGrad: %.6f\n', norm(multirobot_rgrad)); % Changed from Grad
fprintf('  RPE PUDQ: %.6f\n', multirobot_rpe_pudq);
fprintf('  RPE Euclidean: %.6f (Linear: %.6f, Angular: %.6f)\n', multirobot_rpe_eucl, multirobot_rpe_linear, multirobot_rpe_angular);
fprintf('  ATE PUDQ: %.6f\n', multirobot_ate_pudq);
fprintf('  ATE Euclidean: %.6f (Linear: %.6f, Angular: %.6f)\n', multirobot_ate_eucl, multirobot_ate_linear, multirobot_ate_angular);

fprintf('Cartin-sync Results:\n');
fprintf('  SE2 Cost: %.6f\n', cartan_se2_cost);
fprintf('  PUDQ Cost: %.6f\n', cartan_F);
fprintf('  RGrad: %.6f\n', norm(cartan_rgrad)); % Changed from Grad
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

fprintf('\n*** Analysis complete! ***\n');