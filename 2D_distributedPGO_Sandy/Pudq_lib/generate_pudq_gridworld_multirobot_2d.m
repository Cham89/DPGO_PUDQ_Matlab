function graph = generate_pudq_gridworld_multirobot_2d(num_vertices_per_robot, num_robots, sigma_pudq, lc_dist, lc_prob)

graph = {}; 
df = 10; % Wishart df parameter
mu_ij = zeros(3,1); % noise

num_vertices = num_vertices_per_robot * num_robots;

vertices_true = cell(1, num_vertices);
vertices_pudq_true = cell(1, num_vertices);

v_origin = [0.0; 0.0; 0.0];
vertices_true{1} = v_origin;
vertices_pudq_true{1} = pose_to_pudq(v_origin);

for i = 2 : num_vertices
    v_pudq_prev = vertices_pudq_true{i - 1};
    random_delta = pose_to_pudq(randomdelta(0.5));
    v_new = pudq_compose(v_pudq_prev, random_delta);
    vertices_true{i} = pudq_to_pose(v_new);
    vertices_pudq_true{i} = v_new;
end

loop_closures = [];
lc_index = 1;
for i = 1:num_vertices
    for j = 1:num_vertices
        if abs(i - j) > 1
            ran_lc_prob = rand;
            if norm(vertices_true{i} - vertices_true{j}) <= lc_dist && ran_lc_prob < lc_prob
                    loop_closures(:, lc_index) = [i; j];
                    lc_index = lc_index + 1;
            end
        end
    end
end

% Add intra lc for each robot's last and first vertex / and inter from
% robot 1 to each robot
for rr = 1: num_robots
    first_vertex = (rr - 1) * num_vertices_per_robot + 1;
    last_vertex = rr * num_vertices_per_robot;
    second_vertex = (rr - 1) * num_vertices_per_robot + 2;
    third_vertex = (rr - 1) * num_vertices_per_robot + 3;
    loop_closures = lc_search_add(loop_closures, last_vertex, first_vertex);
    loop_closures = lc_search_add(loop_closures, 2, second_vertex);
    loop_closures = lc_search_add(loop_closures, third_vertex, 3);
end

size_lc = size(loop_closures);
num_lc = size_lc(2);
num_edges = num_vertices - 1 + num_lc;

vertices = cell(1, num_vertices);
vertices_pudq = cell(1, num_vertices);

edge_id = 0;
edges = zeros(num_edges, 2);
delta_poses = cell(1, num_edges);
delta_poses_true = cell(1, num_edges);
delta_poses_pudq = cell(1, num_edges);
delta_poses_pudq_true = cell(1, num_edges);

t = cell(1, num_edges);
R = cell(1, num_edges);
information = cell(1, num_edges);
information_pudq = cell(1, num_edges);
tau = cell(1, num_edges);
kappa = cell(1, num_edges);

% Consecutive edges and odom vertices
for i = 1 : num_vertices
    if i == 1
        vertices{i} = zeros(3,1);
        vertices_pudq{i} = pose_to_pudq(vertices{i});
    else
        edge_id = edge_id + 1;
        edges(edge_id, :) = [i-1, i];
        x_i = vertices_pudq_true{i-1};
        x_j = vertices_pudq_true{i};

        Sigma_ij = generate_pudq_cov(sigma_pudq, df);
        Omega_ij = inv(Sigma_ij);
        eta_ij = mvnrnd(mu_ij, Sigma_ij, 1)';
        exp_eta_ij = Exp_1(eta_ij);

        z_ij_true = pudq_compose(pudq_inv(x_i), x_j);
        z_ij = pudq_compose(z_ij_true, exp_eta_ij);

        delta_poses_pudq{edge_id} = z_ij;
        delta_poses_pudq_true{edge_id} = z_ij_true;
        information_pudq{edge_id} = Omega_ij;

        delta_poses_tilde = pudq_to_pose(z_ij);
        delta_poses{edge_id} = delta_poses_tilde;
        delta_poses_true{edge_id} = pudq_to_pose(z_ij_true);

        Omega_ij_eucl = info_pudq_to_eucl(Omega_ij, 2*eta_ij(1)); 
        information{edge_id} = Omega_ij_eucl;
        information_se2{edge_id} = info_pudq_to_se2(Omega_ij); 

        t{edge_id} = delta_poses_tilde(1:2);
        dth_tilde = delta_poses_tilde(3);                      % 要研究此程式轉換
        R{edge_id} = [cos(dth_tilde), -sin(dth_tilde);
                     sin(dth_tilde), cos(dth_tilde)];

        x_i_odom = vertices_pudq{i-1};
        x_j_odom = pudq_compose(x_i_odom, z_ij);
        vertices_pudq{i} = x_j_odom;
        vertices{i} = pudq_to_pose(x_j_odom);

        tau{edge_id} = 2 / trace(inv(Omega_ij_eucl(1:2, 1:2))); 
        kappa{edge_id} = Omega_ij_eucl(3,3);                     
    end 
end

% Non consecutive edges
for i = 1: num_lc
    edge_id = edge_id + 1;
    edges(edge_id, :) = loop_closures(:, i)';

    x_i = vertices_pudq_true{loop_closures(1,i)};
    x_j = vertices_pudq_true{loop_closures(2,i)};

    Sigma_ij = generate_pudq_cov(sigma_pudq, df);
    Omega_ij = inv(Sigma_ij);
    eta_ij = mvnrnd(mu_ij, Sigma_ij, 1)';
    exp_eta_ij = Exp_1(eta_ij);

    z_ij_true = pudq_compose(pudq_inv(x_i), x_j);
    z_ij = pudq_compose(z_ij_true, exp_eta_ij);

    delta_poses_pudq{edge_id} = z_ij;
    delta_poses_pudq_true{edge_id} = z_ij_true;
    information_pudq{edge_id} = Omega_ij;

    delta_poses_tilde = pudq_to_pose(z_ij);
    delta_poses{edge_id} = delta_poses_tilde;
    delta_poses_true{edge_id} = pudq_to_pose(z_ij_true);

    Omega_ij_eucl = info_pudq_to_eucl(Omega_ij, 2*eta_ij(1)); 
    information{edge_id} = Omega_ij_eucl;
    information_se2{edge_id} = info_pudq_to_se2(Omega_ij); 

    t{edge_id} = delta_poses_tilde(1:2);
    dth_tilde = delta_poses_tilde(3);                      
    R{edge_id} = [cos(dth_tilde), -sin(dth_tilde);
                 sin(dth_tilde), cos(dth_tilde)];

    tau{edge_id} = 2 / trace(inv(Omega_ij_eucl(1:2, 1:2)));  
    kappa{edge_id} = Omega_ij_eucl(3,3);                     
end

graph.vertices = vertices;
graph.vertices_pudq = vertices_pudq;
graph.vertices_pudq_true = vertices_pudq_true;
graph.vertices_true = vertices_true;

graph.edges = edges;
graph.delta_poses = delta_poses;
graph.delta_poses_true = delta_poses_true;
graph.delta_poses_pudq = delta_poses_pudq;
graph.delta_poses_pudq_true = delta_poses_pudq_true;
graph.t = t;
graph.R = R;
graph.information = information;
graph.information_pudq = information_pudq;
graph.information_se2 = information_se2;
graph.kappa = kappa;
graph.tau = tau;

graph.Omega = Omega_sparse(graph); 

fprintf('Generated Multi-Robot Square GridWorld with %d robots, %d vertices and %d edges.\n', num_vertices, num_edges);

end

























