function dv_lb = compute_dv_lower_bound(a, d_delta_e, d_delta_i)
% Computes the delta-v lower bound for a maneuver using Table 5.13 from Chernick
%
    % Compute mean motion
    mu = 3.986004418e14; % m^3/s^2
    n = sqrt(mu / (a*1e3)^3);

    % In-plane eccentricity vector change
    ecc_dv_min = norm(a * d_delta_e) * n / 2;

    % Cross-track inclination change
    i_dv_min = norm(a * d_delta_i) * n;

    dv_lb = sqrt(ecc_dv_min^2 + i_dv_min^2);
    %dv_lb = min(ecc_dv_min,i_dv_min);
end
