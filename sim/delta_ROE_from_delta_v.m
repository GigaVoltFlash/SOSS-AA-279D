function delta_roe = delta_ROE_from_delta_v(delta_v_RTN, chief_qns)
    % Computes the change in relative orbital elements (delta_roe)
    % from an impulsive delta-v in RTN frame using quasi-nonsingular elements.
    %
    % Inputs:
    %   delta_v_RTN - 3x1 vector: [dv_R; dv_T; dv_N] (in km/s)
    %   chief_qns   - 6x1 vector: [a, e_x, e_y, i, RAAN, u] (radians)
    %
    % Output:
    %   delta_roe   - 6x1 vector: change in relative orbital elements

    global mu_earth;

    % Unpack QNS elements
    a   = chief_qns(1);
    u_k = deg2rad(chief_qns(6));  % mean argument of latitude

    % Mean motion
    n = sqrt(mu_earth / a^3);

    % Build Gamma_k matrix
    Gamma_k = [  0,        2,           0;
                -2,        0,           0;
             sin(u_k), 2*cos(u_k),      0;
            -cos(u_k), 2*sin(u_k),      0;
                 0,        0,      cos(u_k);
                 0,        0,      sin(u_k)];

    % Scale by 1/(n*a)
    delta_roe = (1 / (n * a)) * Gamma_k * delta_v_RTN;
end