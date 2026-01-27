function log_1_rij = e_ij(z_ij, x_i, x_j)
    r_ij = pudq_compose(pudq_inv(z_ij), pudq_mul(pudq_inv(x_i), x_j));
    log_1_rij = Log_1(r_ij);
end
