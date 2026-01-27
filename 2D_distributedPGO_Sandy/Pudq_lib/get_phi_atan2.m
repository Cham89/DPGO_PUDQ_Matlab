function phi = get_phi_atan2(sin_phi, cos_phi)
   
    %Prevent theta out of range [-pi, pi]
    phi = atan2(sin_phi, cos_phi);
    if phi <= -pi/2
        phi = phi + pi;
    elseif phi > pi/2
        phi = phi - pi;
    end
end
