function exp_x_n = Exp_X_N(X, Y_T)
    num_vertices = length(X)/4;
    exp_x_n = zeros(num_vertices*4, 1);
    
    for i=1:num_vertices
        i_index = 4*(i-1)+1;
        x_i = X(i_index:i_index+3, 1);
        y_t = Y_T(i_index:i_index+3, 1);

        exp_x = Exp_x(x_i, y_t);
        exp_x_n(i_index:i_index+3, 1) = exp_x;
    end
end
