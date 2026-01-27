function G_optimized = reconstructglobalgraph(robot, original_G)

    num_robots = length(robot);
    total_vertices = length(original_G.vertices);

    fprintf('Reconstructing global graph from %d robots (%d total vertices)...\n', num_robots, total_vertices);

    G_optimized = original_G;

    vertices_optimized = cell(1, total_vertices);
    vertices_pudq_optimized = cell(1, total_vertices);

    for r = 1:num_robots
        local_graph = robot(r).local_graph;
        global_mapping = local_graph.vertices_interval;

        fprintf('Robot %d: processing %d local vertices (global indices %d to %d)\n', r, length(local_graph.vertices_pudq), min(global_mapping), max(global_mapping));

        for local_idx = 1:length(local_graph.vertices_pudq)
            if local_idx <= length(global_mapping)
                global_idx = global_mapping(local_idx);
                if global_idx > 0 && global_idx <= total_vertices
                    vertices_optimized{global_idx} = pudq_to_pose(local_graph.vertices_pudq{local_idx});
                    vertices_pudq_optimized{global_idx} = local_graph.vertices_pudq{local_idx};
                end
            end
        end
    end

    G_optimized.vertices = vertices_optimized;
    G_optimized.vertices_pudq = vertices_pudq_optimized;

    fprintf('Global graph reconstruction complete:\n');
end