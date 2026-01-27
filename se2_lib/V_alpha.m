function V = V_alpha(alpha)
    V = [sinc(alpha), -cosc(alpha);
         cosc(alpha), sinc(alpha)];
end
