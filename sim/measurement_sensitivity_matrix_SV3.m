function H = measurement_sensitivity_matrix_SV3(r_chief,v_chief,SV1_oe_qns)
    u_1 = deg2rad(SV1_oe_qns(6));

    J = [...
    1, 0, -cos(u_1), -sin(u_1), 0, 0;
    0, 1,  2*sin(u_1), -2*cos(u_1), 0, 0;
    0, 0,  0, 0, sin(u_1), -cos(u_1)];

    first_entry = J;

    % Build RTN basis vectors
    R_hat = r_chief ./ vecnorm(r_chief, 2, 2);
    h_vec = cross(r_chief, v_chief, 2);
    N_hat = h_vec ./ vecnorm(h_vec, 2, 2);
    T_hat = cross(N_hat, R_hat, 2);

    Q_eci2rtn = [R_hat', T_hat', N_hat']';

    second_entry = Q_eci2rtn' * J;

    H = 0.001*[first_entry;second_entry]; % Convert to kilometers
end

