% === SCRIPT: load_full_analysis_data.m ===

clear, clc, close all;
run('Cartan-Sync\MATLAB\setup.m');
addpath([pwd '/Pudq_lib']);
addpath([pwd '/riem_opt_lib']);
addpath([pwd '/se2_lib']);

% Parameters
num_robots     = 3; % Must match C++ main.cpp
cpp_export_dir = 'cpp_export_full'; % Must match C++ output_dir

fprintf('Loading all data from: %s\n', cpp_export_dir);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 1. LOAD AND REBUILD THE 'graph' STRUCT %%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = struct;

% --- Helper function to load matrix and convert to cell array ---
    function out_cell = load_mat_as_cell(filename)
        mat = readmatrix(filename);
        [N, ~] = size(mat);
        out_cell = cell(N, 1);
        for idx = 1:N
            out_cell{idx} = mat(idx, :)';
        end
    end

% --- Load all graph data ---
G.vertices_pudq       = load_mat_as_cell(fullfile(cpp_export_dir, 'optimized_vertices_pudq.txt'));
G.vertices            = load_mat_as_cell(fullfile(cpp_export_dir, 'optimized_vertices.txt'));
G.vertices_pudq_true  = load_mat_as_cell(fullfile(cpp_export_dir, 'true_vertices_pudq_true.txt'));
G.vertices_true       = load_mat_as_cell(fullfile(cpp_export_dir, 'true_vertices.txt'));

G.delta_poses_pudq_true = load_mat_as_cell(fullfile(cpp_export_dir, 'true_delta_poses_pudq.txt'));
G.delta_poses_true    = load_mat_as_cell(fullfile(cpp_export_dir, 'true_delta_poses.txt'));

% Load edges and convert from 0-indexed to 1-indexed
G.edges = readmatrix(fullfile(cpp_export_dir, 'edges_0_indexed.txt')) + 1;

fprintf('Loaded graph data (N=%d, M=%d)\n', length(G.vertices), length(G.edges));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 2. LOAD AND REBUILD THE 'robot' STRUCT %%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
robot = struct;
for r = 1:num_robots
    robot_idx_c = r - 1; % C++ is 0-indexed
    robot(r).id = r;
    
    cost_file = fullfile(cpp_export_dir, sprintf('robot_%d_cost.txt', robot_idx_c));
    grad_file = fullfile(cpp_export_dir, sprintf('robot_%d_gradnorm.txt', robot_idx_c));
    
    robot(r).rgn_cost = readmatrix(cost_file);
    robot(r).rgn_gradnorm = readmatrix(grad_file);
end
fprintf('Loaded cost/grad history for %d robots.\n', num_robots);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% 3. RUN ALL YOUR MATLAB CODE (UNCHANGED) %%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--- Running MATLAB Analysis ---');

% --- Run ATE/RPE (using the reoriented graph) ---
G_multirobot = reorient_graph(G); % Use your anchor-to-origin function

multirobot_rpe_pudq = RPE_G_pudq(G_multirobot);
[multirobot_rpe_eucl, multirobot_rpe_linear, multirobot_rpe_angular] = RPE_G_eucl(G_multirobot);
multirobot_ate_pudq = ATE_G_pudq(G_multirobot);
[multirobot_ate_eucl, multirobot_ate_linear, multirobot_ate_angular] = ATE_G_eucl(G_multirobot);

fprintf('Multi-Robot Results (from C++ Data):\n');
fprintf('  RPE PUDQ: %.6f\n', multirobot_rpe_pudq);
fprintf('  RPE Euclidean: %.6f (Linear: %.6f, Angular: %.6f)\n', multirobot_rpe_eucl, multirobot_rpe_linear, multirobot_rpe_angular);
fprintf('  ATE PUDQ: %.6f\n', multirobot_ate_pudq);
fprintf('  ATE Euclidean: %.6f (Linear: %.6f, Angular: %.6f)\n', multirobot_ate_eucl, multirobot_ate_linear, multirobot_ate_angular);


% --- Run Plotting ---
figure('Name','Per-Robot Cost (from C++)','NumberTitle','off');
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

figure('Name','Per-Robot Grad Norm (from C++)','NumberTitle','off');
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

disp('Analysis and plotting complete.');