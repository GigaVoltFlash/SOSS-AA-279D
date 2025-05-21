function delta_roe = delta_ROE_from_delta_v(delta_v_RTN, chief_oe)
    % Computes the change in relative orbital elements (delta_roe)
    % from an impulsive delta-v in RTN frame.
    %
    % Inputs:
    %   delta_v_RTN            - 3x1 vector: [dv_R; dv_T; dv_N] (in km/s)
    %   chief_orbital_elements - struct or 1x6 vector: [a, e, i, RAAN,
    %   omega, nu] rad
    %
    % Output:
    %   delta_roe              - 6x1 vector: change in relative orbital elements

    global mu_earth;

    a = chief_oe(1);
    omega = deg2rad(chief_oe(5));
    M = true2mean(deg2rad(chief_oe(6))); % nu (true anomaly) to mean anomaly

    % Mean motion
    n = sqrt(mu_earth / a^3);

    % Mean argument of latitude
    u_k = omega + M;

    % Build transformation matrix Gamma_k
    Gamma_k = [  0,        2,           0;
                -2,        0,           0;
             sin(u_k), 2*cos(u_k),      0;
            -cos(u_k), 2*sin(u_k),      0;
                 0,        0,      cos(u_k);
                 0,        0,      sin(u_k)];

    % Scale by 1/(n*a)
    delta_roe = (1 / (n * a)) * Gamma_k * delta_v_RTN;
end