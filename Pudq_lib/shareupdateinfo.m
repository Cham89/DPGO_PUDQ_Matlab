function robot = shareupdateinfo(share_robot_id, share_subgraph, robotList)

    robot = robotList;
    num_robots = length(robotList);
    % fprintf('[SHARING] Robot %d broadcasting updates to all other robots\n', share_robot_id);

    for recieve_robot_id = 1 : num_robots
        if recieve_robot_id == share_robot_id
            continue;
        end
        recieve_robot_graph = robot(recieve_robot_id).local_graph;
        update_num = 0;

        for L = 1 : length(recieve_robot_graph.lm_vertices)
            landmark_info = recieve_robot_graph.lm_foreign_info{L};
       
            if landmark_info.robot == share_robot_id 
                update_local_idx = landmark_info.vertex;

                if update_local_idx <= length(share_subgraph.vertices)
                    old_pose = recieve_robot_graph.lm_vertices{L};
                    new_pose = share_subgraph.vertices{update_local_idx};
                    new_pose_pudq = share_subgraph.vertices_pudq{update_local_idx};
    
                    recieve_robot_graph.lm_vertices{L} = new_pose;
                    recieve_robot_graph.lm_vertices_pudq{L} = new_pose_pudq;
                    update_num = update_num + 1;
    
                    %fprintf('  Robot %d: Updated landmark %d (R%d-V%d): [%.6f,%.6f,%.6f] -> [%.6f,%.6f,%.6f]\n', ...
                            %recieve_robot_id, L, share_robot_id, update_local_idx, ...
                            %old_pose(1), old_pose(2), old_pose(3), ...
                            %new_pose(1), new_pose(2), new_pose(3));
                else
                    error('Invalid vertex index %d for robot %d (max: %d)', update_local_idx, share_robot_id, length(share_subgraph.vertices));
            
                end
            end
        end

        robot(recieve_robot_id).local_graph = recieve_robot_graph;

        if update_num > 0
            %fprintf('  Robot %d: Applied %d updates from Robot %d\n', recieve_robot_id, update_num, share_robot_id);
        end
    end
    %fprintf('[SHARING] Robot %d update propagation complete\n', share_robot_id);
end






