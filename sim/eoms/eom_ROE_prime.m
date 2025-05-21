function [delta_alpha_qns_prime_dot] = eom_ROE_prime(t, deputy_qns, chief_qns)
    % Computes δα'_qns_dot as a 6x1 vector according to Koenig Eq. 23

    % Inputs: current chief and deputy quasi-nonsing OE (a,e_x,e_yi,RAAN,omega,nu) in
    % deg

    global mu_earth;
    global R_earth;
    global J2;

    % Unpack QNS elements for deputy and chief
    a_d = deputy_qns(1);
    e_x_d = deputy_qns(2);
    e_y_d = deputy_qns(3);
    i_d = deg2rad(deputy_qns(4));

    a_c = chief_qns(1);
    e_x_c = chief_qns(2);
    e_y_c = chief_qns(3);
    i_c = deg2rad(chief_qns(4));

    % Reconstruct magnitudes and angles
    e_d = sqrt(e_x_d^2 + e_y_d^2);
    e_c = sqrt(e_x_c^2 + e_y_c^2);

    omega_d = atan2(e_y_d, e_x_d);
    omega_c = atan2(e_y_c, e_x_c);
    
    % Compute eta
    eta_d = sqrt(1 - e_d^2);
    eta_c = sqrt(1 - e_c^2);

    % Compute kappa
    kappa_d = (3 * J2 * R_earth^2 * sqrt(mu_earth)) / (4 * a_d^(7/2) * eta_d^4);
    kappa_c = (3 * J2 * R_earth^2 * sqrt(mu_earth)) / (4 * a_c^(7/2) * eta_c^4);

    % Trigonometric precomputations
    cos_id = cos(i_d);
    cos_ic = cos(i_c);
    sin_ic = sin(i_c);
    
    cos2_id = cos_id^2;
    cos2_ic = cos_ic^2;
    
    delta_omega = omega_d - omega_c;

    % κ_d contribution
    vec_d = zeros(6,1);
    vec_d(2) = eta_d * (3 * cos2_id - 1) + (5 * cos2_id - 1) - 2 * cos_id * cos_ic;
    vec_d(3) = -e_d * sin(delta_omega) * (5 * cos2_id - 1);
    vec_d(4) =  e_d * cos(delta_omega) * (5 * cos2_id - 1);
    vec_d(6) = -2 * cos_id * sin_ic;

    % -κ_c contribution
    vec_c = zeros(6,1);
    vec_c(2) = (1 + eta_c) * (3 * cos2_ic - 1);
    vec_c(3) = -e_d * sin(delta_omega) * (5 * cos2_ic - 1);
    vec_c(4) =  e_d * cos(delta_omega) * (5 * cos2_ic - 1);
    vec_c(6) = -2 * cos_ic * sin_ic;

    % Final result
    delta_alpha_qns_prime_dot = kappa_d * vec_d - kappa_c * vec_c;
end