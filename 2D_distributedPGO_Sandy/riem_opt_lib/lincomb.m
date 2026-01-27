function v = lincomb(a1, d1, a2, d2)
    if nargin == 2
        v = a1*d1;
    elseif nargin == 4
        v = a1*d1 + a2*d2;
    else
        error('lincomb takes either 2 or 4 inputs.');
    end
end
