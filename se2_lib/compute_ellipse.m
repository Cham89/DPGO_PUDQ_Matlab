function ellipse_xyz = compute_ellipse(phi, lambdas, sigma)
    

%     x = sigma*sqrt(lambdas(1))*cos(phi);
%     y = sigma*sqrt(lambdas(2))*sin(phi);
%     z = sigma*sqrt(lambdas(3))*small_sqrt(1.0 - x^2/(sigma^2*lambdas(1)) - y^2/(sigma^2*lambdas(2)));
    
    y = sigma*sqrt(lambdas(2))*cos(phi);
    z = sigma*sqrt(lambdas(3))*sin(phi);
    x = sigma*sqrt(lambdas(1))*small_sqrt(1.0 - y^2/(sigma^2*lambdas(2)) - z^2/(sigma^2*lambdas(3)));
    
    ellipse_xyz = [x; y; z];
end