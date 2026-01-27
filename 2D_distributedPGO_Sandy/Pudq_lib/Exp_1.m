function exp_1 = Exp_1(x_t)
% Exponential map (tech report eq (5))
    gamma = sinc1(x_t(1));
    q = [cos(x_t(1)), gamma*x_t']';
    exp_1 = pudq_normalize(q);
end