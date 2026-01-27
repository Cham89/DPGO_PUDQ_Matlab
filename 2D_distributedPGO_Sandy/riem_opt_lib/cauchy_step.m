function eta = cauchy_step(H, g, Delta)
    gHg = g'*H*g;
    if gHg > 0
        alpha = min(g'*g/gHg, Delta/norm(g));
    else
        alpha = Delta/norm(g);
    end
    eta = -alpha*g;
end