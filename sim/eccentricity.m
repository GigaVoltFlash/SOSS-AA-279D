% Calculating eccentricity vector given ECI position and velocity
function e = eccentricity(r_IJK,v_IJK)
    mu_earth = 3.986e5; % km^3/s^2
    
    h = cross(r_IJK,v_IJK);
    e = (cross(v_IJK,h)/mu_earth) - (r_IJK/norm(r_IJK));
end