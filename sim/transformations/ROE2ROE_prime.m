function delta_alpha_qns_prime = ROE2ROE_prime(delta_alpha_qns, chief_orbital_elements)
    % Transform delta_alpha_qns to delta_alpha_qns_prime using J_qns(omega_c)
    %
    % Inputs:
    %   delta_alpha_qns         - 6x1 vector in QNS frame
    %   chief_orbital_elements  - struct with field .omega [rad]
    %
    % Output:
    %   delta_alpha_qns_prime   - 6x1 vector in QNS' frame

    % Extract omega_c (argument of perigee) from chief's orbital elements
    omega_c = chief_orbital_elements(5);

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