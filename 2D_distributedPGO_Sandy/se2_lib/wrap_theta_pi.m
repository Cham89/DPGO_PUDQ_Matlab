function theta_wrapped = wrap_theta_pi(theta)
    if theta < -pi
        theta_wrapped = theta + 2*pi;
    elseif theta > pi
        theta_wrapped = theta - 2*pi;
    else
        theta_wrapped = theta;
    end
end
