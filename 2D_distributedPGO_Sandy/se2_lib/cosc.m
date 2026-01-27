function val = cosc(theta)
    order = 12; %n-th order Taylor series approximation
    val = 0;
    for n=1:order
        val = val + ((-1)^(n+1)*theta^(2*n-1))/factorial(2*n);
    end
end
