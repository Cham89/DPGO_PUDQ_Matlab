function se_sync_graph = se_sync_to_graph(graph, xhat)
    %Copy graph structure over
    se_sync_graph = graph;
    
    d = length(xhat.t(:,1)); %Dimension
    n = length(xhat.t(1,:)); %# of vertices
    
    se_sync_vertices = cell(1,n);
    se_sync_vertices_pudq = cell(1,n);
    
    delta_t = zeros(2,1);
    delta_R = eye(2);
    
    %Loop through vertices and copy them back into the graph structure
    for i=1:n
       t_i = xhat.t(:,i);
       R_index = 2*(i-1)+1;
       R_i = xhat.R(:, R_index:R_index+1);
       if i==1
           delta_R = R_i';
           delta_t = -R_i'*t_i;
       end
       
       %Undo the transformation at the origin
       R_i = delta_R*R_i;
       t_i = delta_R*t_i + delta_t;
       
       theta_i = theta_to_R(R_i);
       pose_i = [t_i', theta_i]';
       
       se_sync_vertices{i} = pose_i;
       se_sync_vertices_pudq{i} = pose_to_pudq(pose_i);
    end
    
    se_sync_graph.vertices = se_sync_vertices;
    se_sync_graph.vertices_pudq = se_sync_vertices_pudq;
end