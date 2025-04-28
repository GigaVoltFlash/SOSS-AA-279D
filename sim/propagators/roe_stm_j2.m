function [roe_output] = roe_stm_j2(t, qns_roe_init, SV1_oe_init)
    % Inputs in degrees
    % Input t is the all the time-steps that we want to propagate the ROE states by
    % qns_roe_init is the initial conditions of the deputy in quasi-nonsingular relative orbital elements
    % SV1_oe_init is the initial conditions of the chief in absolute singular orbital elements
    global J2 Re mu_earth
    a_SV1_init = SV1_oe_init(1);
    e_SV1_init = SV1_oe_init(2)
    i_SV1_init = SV1_oe_init(3);
    RAAN_SV1_init = SV1_oe_init(4);
    w_SV1_init = SV1_oe_init(5)
    eta = sqrt(1 - e_SV1_init^2);
    E = 1 + eta;
    kappa = 3/4*J2*Re^2 * sqrt(mu_earth)/(a_SV1_init^(7/2) * eta^4);
    F = 4 + 3*eta;
    G = 1/eta^2;
    P = 3*(cosd(i_SV1_init))^2 - 1;
    Q = 5*(cosd(i_SV1_init))^2 - 1;
    R = cosd(i_SV1_init);
    S = sind(2*i_SV1_init);
    T = (sind(i_SV1_init))^2;
    U = sind(i_SV1_init);
    V = tand(i_SV1_init/2);
    W = (cosd(i_SV1_init/2))^2;
    n = sqrt(mu_earth/a_SV1_init^3);
    omegadot = kappa * Q;

    exi = e_SV1_init*cosd(w_SV1_init);
    eyi = e_SV2_init*cosd(w_SV1_init); 

    roe_output = zeros((n_t, 6));
    roe_output(1, :) = qns_roe_init;

    n_t = length(t);
    for iter = 1:n_t
        tau = t(iter);
        exf = e_SV1_init*cosd(omegadot*tau + w_SV1_init);
        eyf = e_SV2_init*cosd(omegadot*tau + w_SV1_init);
        STM = [1,                            0, 0,                                           0,                                            0,                  0;...
            -(3/2*n + 7/2*kappa*E*P)*tau, 1, kappa*exi*F*G*P*tau,                         kappa*eyi*F*G*P*tau,                          -kappa*F*S*tau,     0;...
            7/2*kappa*eyf*Q*tau,          0, cosd(omega_dot*tau)-4*kappa*exi*eyf*G*Q*tau, -sind(omega_dot*tau)-4*kappa*eyi*eyf*G*Q*tau, 5*kappa*eyf*S*tau,  0;...
            -7/2*kappa*exf*Q*tau,         0, sind(omega_dot*tau)-4*kappa*exi*exf*G*Q*tau, cosd(omega_dot*tau)-4*kappa*eyi*exf*G*Q*tau,  -5*kappa*exf*S*tau, 0;...
            0,                            0, 0,                                           0,                                            1,                  0;...
            7/2*kappa*S*tau,              0, -4*kappa*exi*G*S*tau,                        -4*kappa*eyi*G*S*tau,                         2*kappa*T*tau,      1];
        
        roe_output(iter, :) = STM * qns_roe_init;

end