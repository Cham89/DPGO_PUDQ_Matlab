function theta_1 = wrap_theta(theta_0)
    if theta_0 > pi
        theta_1 = theta_0 - 2*pi;
    elseif theta_0 < -pi
        theta_1 = theta_0 + 2*pi;
    else
        theta_1 = theta_0;
    end
end
