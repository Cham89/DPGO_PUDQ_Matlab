function x = se2_vee(X)
    v_1 = X(1,3);
    v_2 = X(2,3);
    alpha = X(2,1);
    x = [v_1; v_2; alpha];
end