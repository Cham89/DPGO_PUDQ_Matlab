function robot = find_robot_for_vertex(multi_graph, num_robots, global_vertex)
    for r = 1:num_robots
        if ismember(global_vertex, multi_graph(r).vertices_interval)
            robot = r;
            return;
        end
    end
    error('Vertex %d not found in any robot', global_vertex);
end