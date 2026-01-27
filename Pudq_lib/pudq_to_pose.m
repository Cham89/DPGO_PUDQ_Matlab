function p = pudq_to_pose(q)
    theta = 2.0*get_phi_atan2(q(2), q(1));
    R_phi = [q(1), q(2);-q(2), q(1)];
    t = 2.0 * R_phi' * q(3:4);
    p = [t', theta]';
end
