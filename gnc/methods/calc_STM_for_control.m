function STM = calc_STM_for_control(T_opt, tf, SV1_oe_init)
    % Some relevant parameters that are useful
    global J2 R_earth mu_earth
    a = SV1_oe_init(1);
    ecc = SV1_oe_init(2);
    i_SV1_init = SV1_oe_init(3);
    i = deg2rad(i_SV1_init);
    eta = sqrt(1 - ecc^2);
    E = 1 + eta;
    kappa = 3/4*J2*R_earth^2 * sqrt(mu_earth)/(a^(7/2) * eta^4);
    F = 4 + 3*eta;
    P = 3*(cos(i))^2 - 1;
    Q = 5*(cos(i))^2 - 1;
    S = sin(2*i);
    T = (sin(i))^2;
    n = sqrt(mu_earth/a^3);
    omega_dot = kappa * Q;

    tau = tf - T_opt; % delta_u / (n + kappa*(eta*P + Q) ?

    STM = [1,                            0, 0,                  0,                   0,              0;...
           -(3/2*n + 7/2*kappa*E*P)*tau, 1, 0,                  0,                   -kappa*F*S*tau, 0;...
           0,                            0, cos(omega_dot*tau), -sin(omega_dot*tau), 0,              0;...
           0,                            0, sin(omega_dot*tau), cos(omega_dot*tau),  0,              0;...
           0,                            0, 0,                  0,                   1,              0;...
           7/2*kappa*S*tau,              0, 0,                  0,                   2*kappa*T*tau,  1];

end