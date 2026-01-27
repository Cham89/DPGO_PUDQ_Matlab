function q = pudq_compose(x, y)
% (tech report eq(36))
    q = pudq_normalize(pudq_mul(x, y));
end
