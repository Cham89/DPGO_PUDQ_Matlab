function export_graph_for_cpp(G, output_dir)

    fprintf('Exporting complete graph data to %s...\n', output_dir);
    
    if ~exist(output_dir, 'dir')
       mkdir(output_dir);
    end

    N = length(G.vertices);
    M = size(G.edges, 1);
    
    fprintf('Graph has %d vertices and %d edges.\n', N, M);

    % --- Convert Cell Arrays to Matrices ---
    vertices_mat = zeros(N, 3);
    vertices_pudq_mat = zeros(N, 4);
    vertices_true_mat = zeros(N, 3);
    vertices_pudq_true_mat = zeros(N, 4); 
    
    delta_poses_mat = zeros(M, 3);
    delta_poses_true_mat = zeros(M, 3);
    delta_poses_pudq_mat = zeros(M, 4);
    delta_poses_pudq_true_mat = zeros(M, 4);
    t_mat = zeros(M, 2); 
    R_mat = zeros(M * 2, 2); 
    
    info_mat = zeros(M * 3, 3); 
    info_pudq_mat = zeros(M * 3, 3);
    info_se2_mat = zeros(M * 3, 3); 
    
    kappa_vec = zeros(M, 1); 
    tau_vec = zeros(M, 1);

    fprintf('Converting cell arrays to matrices...\n');
    for i = 1:N
        if ~isempty(G.vertices{i})
            vertices_mat(i, :) = G.vertices{i}';
        end
        if ~isempty(G.vertices_pudq{i})
            vertices_pudq_mat(i, :) = G.vertices_pudq{i}';
        end
        if ~isempty(G.vertices_true{i})
            vertices_true_mat(i, :) = G.vertices_true{i}';
        end
        if ~isempty(G.vertices_pudq_true{i}) % <-- ADDED
            vertices_pudq_true_mat(i, :) = G.vertices_pudq_true{i}';
        end
    end
    
    for i = 1:M
        if ~isempty(G.delta_poses{i})
            delta_poses_mat(i, :) = G.delta_poses{i}';
        end
        if ~isempty(G.delta_poses_true{i})
            delta_poses_true_mat(i, :) = G.delta_poses_true{i}';
        end
        if ~isempty(G.delta_poses_pudq{i})
            delta_poses_pudq_mat(i, :) = G.delta_poses_pudq{i}';
        end
         if ~isempty(G.delta_poses_pudq_true{i})
            delta_poses_pudq_true_mat(i, :) = G.delta_poses_pudq_true{i}';
        end
        if ~isempty(G.t{i}) 
            t_mat(i, :) = G.t{i}';
        end
        if ~isempty(G.R{i})
             R_mat((i-1)*2 + 1 : i*2, :) = G.R{i};
        end
        if ~isempty(G.information{i}) 
            info_mat((i-1)*3 + 1 : i*3, :) = G.information{i};
        end
        if ~isempty(G.information_pudq{i})
            info_pudq_mat((i-1)*3 + 1 : i*3, :) = G.information_pudq{i};
        end
        if ~isempty(G.information_se2{i}) 
            info_se2_mat((i-1)*3 + 1 : i*3, :) = G.information_se2{i};
        end
        
        kappa_vec(i) = G.kappa{i};
        tau_vec(i) = G.tau{i};
    end
    
    % --- Write Matrices to Files ---
    
    fprintf('Writing text files...\n');
    
    % Vertices
    writematrix(vertices_mat, [output_dir '/vertices.txt'], 'Delimiter', ' ');
    writematrix(vertices_pudq_mat, [output_dir '/vertices_pudq.txt'], 'Delimiter', ' ');
    writematrix(vertices_true_mat, [output_dir '/vertices_true.txt'], 'Delimiter', ' ');
    writematrix(vertices_pudq_true_mat, [output_dir '/vertices_pudq_true.txt'], 'Delimiter', ' '); % 
    
    % Edges
    writematrix(G.edges, [output_dir '/edges.txt'], 'Delimiter', ' ');
    
    % Measurements
    writematrix(delta_poses_mat, [output_dir '/delta_poses.txt'], 'Delimiter', ' ');
    writematrix(delta_poses_true_mat, [output_dir '/delta_poses_true.txt'], 'Delimiter', ' ');
    writematrix(delta_poses_pudq_mat, [output_dir '/delta_poses_pudq.txt'], 'Delimiter', ' ');
    writematrix(delta_poses_pudq_true_mat, [output_dir '/delta_poses_pudq_true.txt'], 'Delimiter', ' ');
    writematrix(t_mat, [output_dir '/t.txt'], 'Delimiter', ' '); 
    writematrix(R_mat, [output_dir '/R.txt'], 'Delimiter', ' '); 
    
    % Information & Stats
    writematrix(info_mat, [output_dir '/information.txt'], 'Delimiter', ' '); 
    writematrix(info_pudq_mat, [output_dir '/information_pudq.txt'], 'Delimiter', ' ');
    writematrix(info_se2_mat, [output_dir '/information_se2.txt'], 'Delimiter', ' ');
    writematrix(kappa_vec, [output_dir '/kappa.txt'], 'Delimiter', ' '); 
    writematrix(tau_vec, [output_dir '/tau.txt'], 'Delimiter', ' ');

    fprintf('Export complete.\n');
end
