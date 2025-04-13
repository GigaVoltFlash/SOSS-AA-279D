function delta_v_t = sma_maneuver(a,e,r,delta_a)
    mu = 3.986e5; % km^3/s^2

    n = sqrt(mu/a^3);
    delta_v_t = n*r*delta_a/(2*a*sqrt(1-e^2));
end