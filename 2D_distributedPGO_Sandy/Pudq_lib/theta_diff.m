function dt = theta_diff(theta_0, theta_1)
    dt = wrap_theta(angle(exp(theta_0*sqrt(-1))*exp(-theta_1*sqrt(-1))));
end
