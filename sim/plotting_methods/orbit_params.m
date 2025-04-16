% Calculating eccentricity vector given ECI position and velocity
function [a,e,i,RAAN,omega,nu,h,e_vec,epsilon] = orbit_params(r_IJK,v_IJK)
    mu_earth = 3.986e5; % km^3/s^2
    
    [a,e,i,RAAN,omega,nu] = ECI2OE(r_IJK,v_IJK);

    h = cross(r_IJK,v_IJK);
    e_vec = (cross(v_IJK,h)/mu_earth) - (r_IJK/norm(r_IJK));
    epsilon = ((norm(v_IJK)^2)/2) - (mu_earth/norm(r_IJK)); 
end