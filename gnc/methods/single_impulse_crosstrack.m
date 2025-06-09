function [delta_v_vals, delta_v_times] = single_impulse_crosstrack(delta_di, SV1_oe_init, u_SV1_init, t0, tf)
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

    % Change for inclination
    % This formulation doesn't account for J2 effects completely, but 
    % that approximation is considered okay because ix = 0 results in the J2
    % drift in inclination being 0.
    % Plus, accounting for J2 requires an iterative solution.
    % Source: https://slab.sites.stanford.edu/sites/g/files/sbiybj25201/files/media/file/space2016_chernickdamico.pdf 
    delta_v_t3 = [0, 0, n*a*norm(delta_di)];
    t_man_t3 = -100;
    k = 0;
    while (t_man_t3 < t0)
        u_man_t3 = atan2(delta_di(2),delta_di(1)) + k*pi;
        t_man_t3 = (u_man_t3 - deg2rad(u_SV1_init)) /(n + kappa*(eta*P + Q)) + t0;
        k = k+1;
    end

    delta_v_times = t_man_t3;
    delta_v_vals = delta_v_t3;
    if any(delta_v_times > tf)
        fprintf("WARNING! There is a maneuver time that is beyond your control window")
    end
end