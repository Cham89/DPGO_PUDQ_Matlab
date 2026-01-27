function dist = SE2_exp_dist(x1, x2)
    %x and y distances are euclidean
    dx = x1(1) - x2(1);
    dy = x1(2) - x2(2);
    
    %theta distance wraps around
    sin_dist = sin(x1(3) - x2(3));
    cos_dist = cos(x1(3) - x2(3));
    dtheta = atan(sin_dist/cos_dist);
    
    dist = [dx; dy; dtheta];
end