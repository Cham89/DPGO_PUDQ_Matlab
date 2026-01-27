function p = pudq_normalize(q)
    %Make sure sin^2(phi) + cos^2(phi) = 1 -> unit quaternion
    p = q;
    p(1:2) = q(1:2) / norm(q(1:2));
end
