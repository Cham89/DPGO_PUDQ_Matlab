function T = pose_to_SE2(p)
    alpha = wrap_theta_pi(p(3));
    T = [cos(alpha), -sin(alpha), p(1);
         sin(alpha), cos(alpha), p(2);
         0, 0, 1];
end
