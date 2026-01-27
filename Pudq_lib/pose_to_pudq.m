function pudq = pose_to_pudq(pose)
    t = pose(1:2);
    theta = pose(3);
    pudq = zeros(4,1);
    pudq(1) = cos(theta/2);
    pudq(2) = sin(theta/2);
    Ham_product_matrix = [pudq(1), pudq(2); 
                        -pudq(2), pudq(1)];
    pudq(3:4) = (1/2) * Ham_product_matrix * t;
end