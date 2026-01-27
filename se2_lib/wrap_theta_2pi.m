function theta_wrapped = wrap_theta_2pi(theta)
    if theta < 0
        theta_wrapped = theta + 2*pi;
    elseif theta > 2*pi
        theta_wrapped = theta - 2*pi;
    else
        theta_wrapped = theta;
    end
end