function V_inv = V_alpha_inv(alpha)
    V_inv = [xcotx(0.5*alpha), 0.5*alpha;
             -0.5*alpha, xcotx(0.5*alpha)];
end
