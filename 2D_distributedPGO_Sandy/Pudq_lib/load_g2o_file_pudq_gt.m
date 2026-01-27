function graph = load_g2o_file_pudq_gt(g2o_file)
% This function is only for 2D g2o file not 3D

file = fopen (g2o_file, 'r');
edge_id = 0;
vertex_id = 0;
readlines = fgets(file);
line_count = 0;

while ischar(readlines)

    token = strtok(readlines);

    if (strcmp(token, 'VERTEX_SE2'))
        [~, id, x, y, theta] = strread(readlines, '%s %d %f %f %f');
        vertex_id = id + 1;
        vertices_true{vertex_id} = [x, y, theta]';
        vertices_pudq_true{vertex_id} = pose_to_pudq([x, y, theta]');

    elseif (strcmp(token, 'EDGE_SE2'))
        [~, id_v1, id_v2, dx, dy, dtheta, I11, I12, I13, I22, I23, I33] = strread(readlines, '%s %d %d %f %f %f %f %f %f %f %f %f');
        edge_id = edge_id + 1; 
        edges(edge_id, :) = [id_v1 + 1, id_v2 + 1];
        delta_poses{edge_id} = [dx, dy, dtheta]';
        delta_poses_pudq{edge_id} = pose_to_pudq(delta_poses{edge_id});
        t{edge_id} = [dx, dy]';
        R{edge_id} = [cos(dtheta), -sin(dtheta);
                      sin(dtheta), cos(dtheta)];
        measurement_info = [I11, I12, I13;
                            I12, I22, I23;
                            I13, I23, I33];
        information{edge_id} = measurement_info;
        information_pudq{edge_id} = info_eucl_to_pudq(measurement_info, 0);
        information_se2{edge_id} = info_eucl_to_se2(measurement_info, 0);
        tau{edge_id} = 2 / trace(inv(measurement_info(1:2, 1:2)));
        kappa{edge_id} = I33;

    else
        fprintf('No data or finish reading');

    end

    readlines = fgets(file);
    line_count = line_count + 1;

    % Debug if is really reading or not
    if (mod(line_count, 1000) == 0)
        disp(['Reading line' num2str(line_count)])
    end
end

fclose(file);

disp(['Done reading now forming GT vertices'])
num_vertices = max(max(edges));
vertices = cell(1, num_vertices);
vertices_pudq = cell(1, num_vertices);

for i = 1 : num_vertices
    if isempty(vertices_true{i}) || isempty(vertices_pudq_true{i})
        error('Ground truth vertex %d is missing!', i);
    else
        vertices{i} = vertices_true{i};
        vertices_pudq{i} = vertices_pudq_true{i};
    end
end
disp(['Done forming GT vertices'])

graph = {};
graph.vertices = vertices;
graph.vertices_pudq = vertices_pudq;
graph.vertices_true = vertices_true;
graph.vertices_pudq_true = vertices_pudq_true;
graph.edges = edges;
graph.delta_poses = delta_poses;
graph.delta_poses_pudq = delta_poses_pudq;
graph.t = t;
graph.R = R;
graph.information = information;
graph.information_pudq = information_pudq;
graph.information_se2 = information_se2;
graph.kappa = kappa;
graph.tau = tau;
graph.Omega = Omega_sparse(graph);

fprintf('Finished GT Data Loading')

end