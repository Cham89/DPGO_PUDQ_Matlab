function alpha = SE2_get_alpha(T)
    alpha = atan2(T(2,1), T(1,1));

    %Wrap alpha from -pi to pi
    alpha = wrap_theta_pi(alpha);
end
