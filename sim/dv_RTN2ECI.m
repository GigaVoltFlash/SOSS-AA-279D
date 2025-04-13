function delta_v_ECI = dv_RTN2ECI(r_ECI, v_ECI, delta_v_RTN)
    % Normalize R = radial direction (along position)
    R_hat = r_ECI / norm(r_ECI);

    % N = orbit normal direction = r × v
    N = cross(r_ECI, v_ECI);
    N_hat = N / norm(N);

    T_hat = cross(N_hat, R_hat);

    % Construct rotation matrix from RTN to ECI
    Q_RTN2ECI = [R_hat, T_hat, N_hat];

    delta_v_ECI = Q_RTN2ECI * delta_v_RTN;
end