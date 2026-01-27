function deltapose = randomdelta(straight_prob)

    random_prob = rand;
    turn_prob = straight_prob + (1-straight_prob)/2;
    if random_prob <= straight_prob
        deltapose = [1.0; 0.0; 0.0];
    elseif random_prob <= turn_prob
        deltapose = [0.0; 1.0; pi/2];
    else
        deltapose = [0.0; -1.0; -pi/2];
    end
end
