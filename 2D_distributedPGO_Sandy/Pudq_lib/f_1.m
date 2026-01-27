function y = f_1(x)
    if x == 0
        y = 0;
    else
        y = csc(x)^2*(sin(x) - x*cos(x));
    end
end
