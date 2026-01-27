function X = se2_hat(x)
    X = [0, -x(3), x(1); x(3), 0, x(2); 0, 0, 0];
end