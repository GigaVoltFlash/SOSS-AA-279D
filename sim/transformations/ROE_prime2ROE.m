function delta_alpha_qns = ROE_prime2ROE(delta_alpha_qns_prime, chief_orbital_elements)
    % Inverse of ROE2ROE_prime: converts from QNS' frame back to QNS frame
    %
    % Inputs:
    %   delta_alpha_qns_prime   - 6x1 vector in QNS' frame
    %   chief_orbital_elements  - vector with 5th element = omega_c [rad]
    %
    % Output:
    %   delta_alpha_qns         - 6x1 vector in QNS frame

    % Extract omega_c (argument of perigee) from chief's orbital elements
    omega_c = chief_orbital_elements(5);

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
