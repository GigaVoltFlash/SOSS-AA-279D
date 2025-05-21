function delta_alpha_qns_prime = ROE2ROE_prime(delta_alpha_qns, chief_qns)
    % Transforms delta_alpha_qns to delta_alpha_qns_prime using J_qns(omega_c)
    %
    % Inputs:
    %   delta_alpha_qns  - 6x1 vector in QNS frame
    %   chief_qns        - 6x1 vector: [a, e_x, e_y, i, RAAN, u]
    %
    % Output:
    %   delta_alpha_qns_prime - 6x1 vector in QNS' frame

    % Extract e_x and e_y from chief's QNS elements
    e_x_c = chief_qns(2);
    e_y_c = chief_qns(3);

    % Reconstruct omega_c from e_x and e_y
    omega_c = atan2(e_y_c, e_x_c);

    % Build J_qns(omega_c)
    I2 = eye(2);
    Z2 = zeros(2);

    cos_omega = cos(omega_c);
    sin_omega = sin(omega_c);

    R = [cos_omega, sin_omega;
        -sin_omega, cos_omega];

    J = [I2, Z2, Z2;
         Z2, R,  Z2;
         Z2, Z2, I2];

    % Apply transformation
    delta_alpha_qns_prime = J * delta_alpha_qns;
end