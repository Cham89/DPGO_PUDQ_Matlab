function theta = dist_theta(theta0, theta1)
    cos_dist = cos(theta0) - cos(theta1);
    sin_dist = sin(theta0) - sin(theta1);
    theta = wrap_theta(atan(sin_dist/cos_dist));
end
