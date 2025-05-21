function delta_alpha_qns = ROE_prime2ROE(delta_alpha_qns_prime, chief_qns)
    % Inverse of ROE2ROE_prime: converts from QNS' frame back to QNS frame
    % using chief orbital elements in quasi-nonsingular form
    %
    % Inputs:
    %   delta_alpha_qns_prime - 6x1 vector in QNS' frame
    %   chief_qns             - 6x1 vector: [a, e_x, e_y, i, RAAN, u]
    %
    % Output:
    %   delta_alpha_qns       - 6x1 vector in QNS frame

    % Extract e_x and e_y from chief's QNS elements
    e_x_c = chief_qns(2);
    e_y_c = chief_qns(3);

    % Reconstruct omega_c from e_x and e_y
    omega_c = atan2(e_y_c, e_x_c);

    % Build inverse J_qns(-omega_c)
    I2 = eye(2);
    Z2 = zeros(2);
    cos_omega = cos(-omega_c);
    sin_omega = sin(-omega_c);

    R_inv = [cos_omega, sin_omega;
            -sin_omega, cos_omega];

    J_inv = [I2, Z2, Z2;
             Z2, R_inv, Z2;
             Z2, Z2, I2];

    % Apply inverse transformation
    delta_alpha_qns = J_inv * delta_alpha_qns_prime;
end
