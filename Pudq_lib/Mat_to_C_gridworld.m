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
num_robots             = 3;
num_vertices_per_robot = 10;
lc_dist                = 3.0;
lc_prob                = 0.02;
sigma_pudq             = 1e-3;
lambda                 = 1;
alpha_init             = 1e-1;
grad_tol               = 1e-2;
lambda_reg             = 1e-3;
export_directory       = 'C:\temp\matlab_export_c';

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
fprintf('Chordal Initialization complete.\n');

export_graph_for_cpp(G, export_directory);
fprintf('MATLAB processing finished. Data is ready for C++.\n');