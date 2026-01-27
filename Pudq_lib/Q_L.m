function Q = Q_L(x)
% (tech report eq(36))
    Q = [
        x(1), -x(2), 0.0,  0.0;
        x(2), x(1),  0.0,  0.0;
        x(3), x(4),  x(1), -x(2);
        x(4), -x(3), x(2), x(1)
    ];
end
