function T_inv = SE2_inv(T)
    R = T(1:2,1:2);
    t = T(1:2, 3);
    T_inv = [R', -R'*t; 0, 0, 1];
end


