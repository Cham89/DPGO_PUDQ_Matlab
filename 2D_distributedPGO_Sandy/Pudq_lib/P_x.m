function P = P_x(x)
    x_r = x(1:2);
    P = [eye(2) - x_r*x_r', zeros(2,2); zeros(2,2), eye(2)];
end