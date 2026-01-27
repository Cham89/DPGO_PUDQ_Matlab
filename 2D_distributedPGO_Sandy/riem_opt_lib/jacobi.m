function L = jacobi(H)
    [N, ~] = size(H);
    D = spdiags(H, 0);
    Dinv = sqrt(1 ./ D);
    L = full(spdiags(Dinv, 0, N, N));
end