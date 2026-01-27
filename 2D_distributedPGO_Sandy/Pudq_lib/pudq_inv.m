function q_inv = pudq_inv(q)
% Pudq inverse (tech report eq(36) above)
    q_inv = [q(1); -q(2); -q(3); -q(4)];
end