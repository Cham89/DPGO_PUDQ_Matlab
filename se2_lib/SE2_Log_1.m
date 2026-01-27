function X = SE2_Log_1(T)
    t = T(1:2, 3);
    cos_alpha = T(1,1);
    sin_alpha = T(2,1);
    alpha = atan2(sin_alpha, cos_alpha);
    
    v = V_alpha_inv(alpha) * t;
    
    X = [[0, -alpha;alpha, 0], v;0, 0, 1];
end
