function val = small_sqrt(x)
    if abs(x) <= 1e-6
        val = 0.0;
    else
        val = sqrt(x) ;
    end
end