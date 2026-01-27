function local_idx = find_local_vertex_index (multi_graph, global_index)
    vertices_interval = multi_graph.vertices_interval;
    global_idx = find(vertices_interval == global_index);
    
    if isempty(global_idx)
        error('Global vertex %d not found in robot mapping', global_vertex);
    end

    local_idx = global_idx(1); 
end
