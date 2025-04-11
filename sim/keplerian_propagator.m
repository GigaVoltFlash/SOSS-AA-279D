function [r_ECI, v_ECI] = keplerian_propagator(a,e,i,RAAN,omega, t_vals, M0)
    % Inputs are always deg
    global mu_earth
    n = sqrt(mu_earth/a^3);
    r_ECI = zeros(length(t_vals), 3);
    v_ECI = zeros(length(t_vals), 3);
    M_vals = deg2rad(M0) + n*(t_vals);
    for j = 1:length(t_vals)
        E = Newton_Raphson(M_vals(j), e, 1e-5);
        nu = rad2deg(2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e)));
        [r_ECI(j, :), v_ECI(j, :)] = OE2ECI(a, e, i, RAAN, omega, nu);
    end
end


