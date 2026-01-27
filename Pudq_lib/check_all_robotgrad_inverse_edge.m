function [robotList_out] = check_all_robotgrad_inverse_edge(robotList, grad_tol)
% This function is for checking if all the robot gradient converge
    
    robotList_out = robotList;
    num_robots = length(robotList);

    for r = 1:num_robots
        
        local_graph_r = robotList_out(r).local_graph;
        anchor_first_vertex = isfield(local_graph_r,'anchor_first') && local_graph_r.anchor_first && isfield(local_graph_r,'robot_id') && local_graph_r.robot_id == 1;
        num_vertices = length(local_graph_r.vertices);

        if anchor_first_vertex
            movable_indices = 5 : 4 * num_vertices;
        else
            movable_indices = 1 : 4 * num_vertices;
        end 
            
        
        [rgrad, ~, ~] = rgn_gradhess_multi_J_forinverse(local_graph_r);
        norm_grad = norm(rgrad(movable_indices));

        if norm_grad < grad_tol
            robotList_out(r).converge = true;
            fprintf('Robot %d converged (grad = %.3g).\n', r, norm_grad);
        else
            robotList_out(r).converge = false;
            fprintf('Robot %d not converged (grad = %.3g).\n', r, norm_grad);
        end

    end

end