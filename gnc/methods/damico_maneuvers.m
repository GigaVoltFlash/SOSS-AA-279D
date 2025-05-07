function [delta_v_vals, delta_v_times] = damico_maneuvers(delta_da, delta_de, delta_di, SV1_oe_init, u_SV1_init, t0, tf)
    % If you want a change in ecc and delta a, then you can do the two-impulse burn.
    global J2 R_earth mu_earth
    a = SV1_oe_init(1);
    ecc = SV1_oe_init(2);
    i_SV1_init = SV1_oe_init(3);
    i = deg2rad(i_SV1_init);
    eta = sqrt(1 - ecc^2);
    kappa = 3/4*J2*R_earth^2 * sqrt(mu_earth)/(a^(7/2) * eta^4);
    P = 3*(cos(i))^2 - 1;
    Q = 5*(cos(i))^2 - 1;
    n = sqrt(mu_earth/a^3);

    % Two-impulse burn for the delta a and eccentricity change
    delta_v_t1 = [0, n*a/4 * ((delta_da) + norm(delta_de)), 0];
    delta_v_t2 = [0, n*a/4 * ((delta_da) - norm(delta_de)), 0];
    t_man_t1 = -100;
    t_man_t2 = -100;
    k = 0;
    while (t_man_t1 < 0) || (t_man_t2 < 0)
        u_man_t1 = atan2(delta_de(2),delta_de(1)) + k*pi;
        u_man_t2 = u_man_t1 + (k+1)*pi;
        t_man_t1 = (u_man_t1 - deg2rad(u_SV1_init)) /(n + kappa*(eta*P + Q)) + t0;
        t_man_t2 = (u_man_t2 - deg2rad(u_SV1_init)) /(n + kappa*(eta*P + Q)) + t0;
        k = k+1;
    end

    
    % Change for inclination
    delta_v_t3 = [0, 0, n*a*norm(delta_di)];
    t_man_t3 = -100;
    k = 0;
    while (t_man_t3 < 0)
        u_man_t3 = atan2(delta_di(2),delta_di(1)) + k*pi;
        t_man_t3 = (u_man_t3 - deg2rad(u_SV1_init)) /(n + kappa*(eta*P + Q)) + t0;
        k = k+1;
    end

    delta_v_vals = [delta_v_t1; delta_v_t2; delta_v_t3];
    delta_v_times = [t_man_t1, t_man_t2, t_man_t3];

    if any(delta_v_times) > tf
        fprintf("WARNING! There is a maneuver time that is beyond your control window")
    end
end