clc;
fprintf('*** RESUMING SIMULATION ***\n');

addpath([pwd '/Pudq_lib']);
addpath([pwd '/riem_opt_lib']);
addpath([pwd '/se2_lib']);


load('C:\Users\hchen804\OneDrive - Georgia Institute of Technology\Desktop\2D_distributedPGO_Sandy\result_ws\half_Grid1000_1_needtocontinue_2.07_2.01_2robots.mat');
fprintf('Loaded state from Iteration %d\n', iter);

while true
    iter = iter + 1;
    [global_time, r] = min(robot_time);

    [local_graph, F_now, norm_grad_F] = Optimization_RLM_inverse_edge(robot(r).local_graph, mu_reg, grad_tol); 

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

fprintf('\n*** SIMULATION CONVERGED! ***\n');


save('converged_robot_state.mat', 'robot');
fprintf('Saved the converged robot data to "converged_robot_state.mat"\n');

fprintf('\nTo run the final analysis, you must now run the')
fprintf('\nplotting/results section of your main script after')
fprintf('\nloading "converged_robot_state.mat" and your original "G" and "G_Cartan" graphs.\n');