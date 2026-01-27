function graph = generate_pudq_squareworld_multirobot_2d(num_robots, sigma_pudq, side_len , gap, inter_lc_prob)

graph = {}; 
df = 10; % Wishart df parameter
mu_ij = zeros(3,1); % noise

num_vertices_per_robot = 4;
num_vertices = num_robots * num_vertices_per_robot;

global_vertices_idx = @(r,k)(r-1)*num_vertices_per_robot + k;

vertices_true = cell(1, num_vertices);
vertices_pudq_true = cell(1, num_vertices);

for r = 1:num_robots
    base_x = (r-1) * (side_len + gap);
    base_y = 0;

    V1 = [base_x; base_y; 0];
    V2 = [base_x + side_len; base_y; 0];
    V3 = [base_x + side_len; base_y + side_len; pi/2];
    V4 = [base_x; base_y + side_len; pi];

    vertices_true{global_vertices_idx(r,1)} = V1;
    vertices_true{global_vertices_idx(r,2)} = V2;
    vertices_true{global_vertices_idx(r,3)} = V3;
    vertices_true{global_vertices_idx(r,4)} = V4;
    
    vertices_pudq_true{global_vertices_idx(r,1)} = pose_to_pudq(V1);
    vertices_pudq_true{global_vertices_idx(r,2)} = pose_to_pudq(V2);
    vertices_pudq_true{global_vertices_idx(r,3)} = pose_to_pudq(V3);
    vertices_pudq_true{global_vertices_idx(r,4)} = pose_to_pudq(V4);
end

loop_closures = [];
lc_index = 1;

% Intra edges ( 4 -> 1 + 2 -> 4 ) 在考慮要不要加1到3
for r = 1: num_robots
    V1 = global_vertices_idx(r,1);
    V2 = global_vertices_idx(r,2);
    V3 = global_vertices_idx(r,3);
    V4 = global_vertices_idx(r,4);
    loop_closures(:, lc_index) = [V4; V1];
    lc_index = lc_index + 1;
    loop_closures(:, lc_index) = [V2; V4];
    lc_index = lc_index + 1;
end

% Inter edges (guarantee at least one inter edges for each robot)
for r1 = 1 : num_robots
    for r2 = 1 : num_robots
        if r1 ~= r2
            k1 = randi(4);
            k2 = randi(4);
            while r2 == r1 + 1 && k1 == 4 && k2 ==1
                k2 = randi(4);
            end
            i = global_vertices_idx(r1, k1);
            j = global_vertices_idx(r2, k2);
            loop_closures(:, lc_index) = [i; j];
            lc_index = lc_index + 1;
                end
    end
end

% Inter edges (extra parts)
for r1 = 1 : num_robots
    for r2 = 1: num_robots
        if r1 ~= r2
            for k1 = 1:4
                for k2 = 1:4
                    if r2 == r1 + 1 && k1 == 4 && k2 ==1
                        continue;
                    end
                    if rand < inter_lc_prob
                        i = global_vertices_idx(r1, k1);
                        j = global_vertices_idx(r2, k2);
                        exist = false;
                        for existlc = 1 : size(loop_closures, 2)
                            if loop_closures(1, existlc) == i && loop_closures(2, existlc) == j
                                exist = true;
                                break;
                            end
                        end
                        if ~exist
                            loop_closures(:, lc_index) = [i; j];
                            lc_index = lc_index + 1;
                        end
                    end
                end
            end
        end
    end
end

size_lc = size(loop_closures);
num_lc = size_lc(2);
num_edges = num_vertices-1 + num_lc;

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
        dth_tilde = delta_poses_tilde(3);                      
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



























