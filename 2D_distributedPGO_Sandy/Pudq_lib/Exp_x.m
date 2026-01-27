function exp_x = Exp_x(x, y_t)
    %Project y_t into the tangent space at x
    y_t = P_x(x) * y_t;
    x_inv_y_t = pudq_mul(pudq_inv(x), y_t);
    
    exp_1 = Exp_1(x_inv_y_t(2:4));
    exp_x = pudq_compose(x, exp_1);
end
