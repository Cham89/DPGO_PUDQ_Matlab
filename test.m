clc;
fprintf('*** RUNNING DIAGNOSTIC TEST ***\n');

% --- SETUP ---
addpath([pwd '/Pudq_lib']);
addpath([pwd '/riem_opt_lib']);
addpath([pwd '/se2_lib']);

lambda = 1; 


load('debug_state.mat');
fprintf('Loaded state from Iteration %d\n', iter);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_reg = 1e-4; 
% lambda_reg = 1e-1;  
% lambda_reg = 1.0;   

fprintf('TESTING WITH NEW lambda_reg = %.1e\n\n', lambda_reg);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


start_iter = iter; 

while true
    iter = iter + 1;
    [global_time, r] = min(robot_time);

    [local_graph, F_now, norm_grad_F] = Optimization_RLM_inverse_edge(robot(r).local_graph, lambda_reg, grad_tol);

    robot(r).local_graph     = local_graph;
    robot(r).rgn_cost(end+1) = F_now;
    robot(r).rgn_gradnorm(end+1) = norm_grad_F;

    robot = shareupdateinfo(r, local_graph, robot);
    robot_time(r) = global_time + exprnd(1/lambda);

    fprintf('Iter=%3d -> Robot %d updated  (F=%.3g, grad=%.3g)\n', iter, r, F_now, norm_grad_F);

    robot = check_all_robotgrad_inverse_edge(robot, grad_tol);

    if all([robot.converge])
        fprintf('All %d robots converged -> terminating at iter=%d.\n', numel(robot), iter);
        break;
    end
end

fprintf('\n*** TEST COMPLETE ***\n');