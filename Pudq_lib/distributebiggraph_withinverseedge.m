function multi_graph = distributebiggraph_withinverseedge(graph, num_robots)
% Distribute big graph into multi robots' subgraph
% both robot have info for inter edges 
multi_graph = struct;

vertices_total = length(graph.vertices);
vertices_per_robot = floor(vertices_total/num_robots);

for r = 1 : num_robots
    start_idx = (r-1)* vertices_per_robot + 1;
    if r == num_robots
        end_idx = vertices_total;
        vertices_num = end_idx - start_idx + 1;
    else
        end_idx = r* vertices_per_robot;
        vertices_num = vertices_per_robot;
    end

    multi_graph(r).vertices_true = cell(1, vertices_num);
    multi_graph(r).vertices_pudq_true = cell(1, vertices_num);
    multi_graph(r).vertices = cell(1, vertices_num);
    multi_graph(r).vertices_pudq = cell(1, vertices_num);

    for i = 1 : vertices_num
        global_vertices = (r-1)* vertices_per_robot + i; % 可能會出錯的地方
        multi_graph(r).vertices_true{i} = graph.vertices_true{global_vertices};
        multi_graph(r).vertices_pudq_true{i} = graph.vertices_pudq_true{global_vertices};
        multi_graph(r).vertices{i} = graph.vertices{global_vertices};
        multi_graph(r).vertices_pudq{i} = graph.vertices_pudq{global_vertices};
    end
    
    multi_graph(r).vertices_interval = start_idx : end_idx;

    multi_graph(r).intra_edges = [];
    multi_graph(r).intra_dp = {};
    multi_graph(r).intra_dp_true = {};
    multi_graph(r).intra_dp_pudq = {};
    multi_graph(r).intra_dp_pudq_true = {};
    multi_graph(r).intra_t = {};
    multi_graph(r).intra_R = {};
    multi_graph(r).intra_info = {};
    multi_graph(r).intra_info_pudq = {};
    multi_graph(r).intra_info_se2 = {};
    multi_graph(r).intra_tau = {};
    multi_graph(r).intra_kappa = {};

    multi_graph(r).inter_edges = [];
    multi_graph(r).inter_dp = {};
    multi_graph(r).inter_dp_true = {};
    multi_graph(r).inter_dp_pudq = {};
    multi_graph(r).inter_dp_pudq_true = {};
    multi_graph(r).inter_t = {};
    multi_graph(r).inter_R = {};
    multi_graph(r).inter_info = {};
    multi_graph(r).inter_info_pudq = {};
    multi_graph(r).inter_info_se2 = {};
    multi_graph(r).inter_tau = {};
    multi_graph(r).inter_kappa = {};

    multi_graph(r).lm_vertices = {};
    multi_graph(r).lm_vertices_pudq = {};
    multi_graph(r).lm_vertices_true = {};
    multi_graph(r).lm_vertices_pudq_true = {};
    multi_graph(r).lm_foreign_info = {};

end

for e = 1 : size(graph.edges, 1)
    v_i = graph.edges(e, 1);
    v_j = graph.edges(e, 2);

    r_i = find_robot_for_vertex(multi_graph, num_robots, v_i);
    r_j = find_robot_for_vertex(multi_graph, num_robots, v_j);

    % Intra edges
    if r_i == r_j
        r = r_i;
        v_i_local = find_local_vertex_index (multi_graph(r), v_i);
        v_j_local = find_local_vertex_index (multi_graph(r), v_j);

        multi_graph(r).intra_edges = [multi_graph(r).intra_edges, [v_i_local; v_j_local]];
        multi_graph(r).intra_dp{end+1} = graph.delta_poses{e};
        multi_graph(r).intra_dp_true{end+1} = graph.delta_poses_true{e};
        multi_graph(r).intra_dp_pudq{end+1} = graph.delta_poses_pudq{e};
        multi_graph(r).intra_dp_pudq_true{end+1} = graph.delta_poses_pudq_true{e};
        multi_graph(r).intra_t{end+1} = graph.t{e};
        multi_graph(r).intra_R{end+1} = graph.R{e};
        multi_graph(r).intra_info{end+1} = graph.information{e};
        multi_graph(r).intra_info_pudq{end+1} = graph.information_pudq{e};
        multi_graph(r).intra_info_se2{end+1} = graph.information_se2{e};
        multi_graph(r).intra_tau{end+1} = graph.tau{e};
        multi_graph(r).intra_kappa{end+1} = graph.kappa{e};

    % Inter edges
    else
        v_i_local = find_local_vertex_index (multi_graph(r_i), v_i);
        v_j_local = find_local_vertex_index (multi_graph(r_j), v_j);

        [Lm_idx_1, Exists_1] = findOrAddLandmark(multi_graph(r_i), r_j, v_j_local);
        if ~Exists_1
            multi_graph(r_i).lm_vertices{Lm_idx_1} = graph.vertices{v_j};
            multi_graph(r_i).lm_vertices_pudq{Lm_idx_1} = graph.vertices_pudq{v_j};
            multi_graph(r_i).lm_vertices_true{Lm_idx_1} = graph.vertices_true{v_j};
            multi_graph(r_i).lm_vertices_pudq_true{Lm_idx_1} = graph.vertices_pudq_true{v_j};
            multi_graph(r_i).lm_foreign_info{Lm_idx_1} = struct('robot', r_j, 'vertex', v_j_local, 'global_vertex', v_j);
        end

        [Lm_idx_2, Exists_2] = findOrAddLandmark(multi_graph(r_j), r_i, v_i_local);
        if ~Exists_2
            multi_graph(r_j).lm_vertices{Lm_idx_2} = graph.vertices{v_i};
            multi_graph(r_j).lm_vertices_pudq{Lm_idx_2} = graph.vertices_pudq{v_i};
            multi_graph(r_j).lm_vertices_true{Lm_idx_2} = graph.vertices_true{v_i};
            multi_graph(r_j).lm_vertices_pudq_true{Lm_idx_2} = graph.vertices_pudq_true{v_i};
            multi_graph(r_j).lm_foreign_info{Lm_idx_2} = struct('robot', r_i, 'vertex', v_i_local, 'global_vertex', v_i);
        end

        multi_graph(r_i).inter_edges = [multi_graph(r_i).inter_edges, [v_i_local; Lm_idx_1; false]];
        multi_graph(r_i).inter_dp{end+1} = graph.delta_poses{e};
        multi_graph(r_i).inter_dp_true{end+1} = graph.delta_poses_true{e};
        multi_graph(r_i).inter_dp_pudq{end+1} = graph.delta_poses_pudq{e};
        multi_graph(r_i).inter_dp_pudq_true{end+1} = graph.delta_poses_pudq_true{e};
        multi_graph(r_i).inter_t{end+1} = graph.t{e};
        multi_graph(r_i).inter_R{end+1} = graph.R{e};
        multi_graph(r_i).inter_info{end+1} = graph.information{e};
        multi_graph(r_i).inter_info_pudq{end+1} = graph.information_pudq{e};
        multi_graph(r_i).inter_info_se2{end+1} = graph.information_se2{e};
        multi_graph(r_i).inter_tau{end+1} = graph.tau{e};
        multi_graph(r_i).inter_kappa{end+1} = graph.kappa{e};

        % 現在使用的是原本的資訊並沒有倒過來 no inverse
        % z_ij_true = graph.delta_poses_pudq_true{e};
        % z_ij = graph.delta_poses_pudq{e};
        % z_ij_true_inv = pudq_inv(z_ij_true);
        % z_ij_inv = pudq_inv(z_ij);
        % pose_inv = pudq_to_pose(z_ij_inv);
        % t_inv = pose_inv(1:2);
        % theta_inv = pose_inv(3);
        % R_inv = [cos(theta_inv), -sin(theta_inv);
        %          sin(theta_inv), cos(theta_inv)];

        multi_graph(r_j).inter_edges = [multi_graph(r_j).inter_edges, [Lm_idx_2; v_j_local; true]];
        multi_graph(r_j).inter_dp{end+1} = graph.delta_poses{e};
        multi_graph(r_j).inter_dp_true{end+1} = graph.delta_poses_true{e};
        multi_graph(r_j).inter_dp_pudq{end+1} = graph.delta_poses_pudq{e};
        multi_graph(r_j).inter_dp_pudq_true{end+1} = graph.delta_poses_pudq_true{e};
        multi_graph(r_j).inter_t{end+1} = graph.t{e};
        multi_graph(r_j).inter_R{end+1} = graph.R{e};
        multi_graph(r_j).inter_info{end+1} = graph.information{e};
        multi_graph(r_j).inter_info_pudq{end+1} = graph.information_pudq{e};
        multi_graph(r_j).inter_info_se2{end+1} = graph.information_se2{e};
        multi_graph(r_j).inter_tau{end+1} = graph.tau{e};
        multi_graph(r_j).inter_kappa{end+1} = graph.kappa{e};
    end
end

for r = 1:num_robots
    multi_graph(r).Omega_ij = Omega_sparse_Multi(multi_graph(r));
end

fprintf('Successfully distributed gridworld:\n');
end



