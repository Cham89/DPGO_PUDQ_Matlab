function p = SE2_to_pose(T)
    t = T(1:2, 3);
    alpha = SE2_get_alpha(T);
    p = [t', alpha]';
end
