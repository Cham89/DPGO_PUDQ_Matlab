function q = SE2_to_pudq(T)
    t = T(1:2,3);
    theta = wrap_theta(SE2_get_alpha(T));

    q = zeros(4,1);
    q(1) = cos(0.5*theta);
    q(2) = sin(0.5*theta);
    Q_r = [q(1), q(2); -q(2), q(1)];
    q(3:4) = 0.5 * Q_r * t;

    q = sign(q(1))*q;
end
