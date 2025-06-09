function [delta_v_vals, delta_v_times] = radial_impulses_inplane(delta_da, delta_de, delta_dlambda, SV1_oe_init, u_SV1_init, t0, tf)
    % If you want a change in ecc without affecting a, then you can do the two-impulse radial burn.
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
    u0 = deg2rad(u_SV1_init); % Converting initial mean argument of latitude to rad.

    % Two-impulse burn for the delta a and eccentricity change
    delta_v_r1 = [n*a/2 * (-(delta_dlambda)/2 + norm(delta_de)), 0, 0];
    delta_v_r2 = [n*a/2 * (-(delta_dlambda)/2 - norm(delta_de)), 0, 0];
    t_man_t1 = -100;
    t_man_t2 = -100;
    k = 0;
    while (t_man_t1 < t0) || (t_man_t2 < t0)
        u1 = atan2(delta_de(1),delta_de(2)) + k*pi;
        u2 = u1 + (k+2)*pi;
        t_man_t1 = convert_u_to_t(u1);
        t_man_t2 = convert_u_to_t(u2);
        k = k+1;
    end

    delta_v_vals = [delta_v_r1; delta_v_r2];
    delta_v_times = [t_man_t1, t_man_t2];

    if any(delta_v_times > tf)
        fprintf("WARNING! There is a maneuver time that is beyond your control window")
    end
    
    function t = convert_u_to_t(u)
        t = (u - u0) /(n + kappa*(eta*P + Q)) + t0;
    end
end