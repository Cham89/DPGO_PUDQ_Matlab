function y = sinc1(x)
    if x == 0
        y = 1;
    else
        y = sin(x)/x;
    end
end
