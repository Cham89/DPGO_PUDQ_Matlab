function x_t = Log_1(r)
    phi = get_phi_atan2(r(2), r(1));
    gamma = sinc1(phi);
    x_t = [r(2)/gamma; r(3)/gamma; r(4)/gamma];
end