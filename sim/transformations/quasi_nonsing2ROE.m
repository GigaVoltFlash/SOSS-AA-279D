function [d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y] = quasi_nonsing2ROE(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o, ...
    a_t, e_x_t, e_y_t, i_t, RAAN_t, u_t)

    a_o_m = a_o * 1e3; % km --> m
    a_t_m = a_t * 1e3; % km --> m

    d_a = a_t_m - a_o_m; % m

    d_e_x = (e_x_t - e_x_o) * a_o_m; % m
    d_e_y = (e_y_t - e_y_o) * a_o_m; % m

    delta_i = mod(i_t - i_o + 180, 360) - 180; % Wrap into [-180,180]
    delta_RAAN = mod(RAAN_t - RAAN_o + 180, 360) - 180;
    delta_u = mod(u_t - u_o + 180, 360) - 180;

    d_i_x = deg2rad(delta_i) * a_o_m; % m
    d_i_y = deg2rad(delta_RAAN) * a_o_m * sind(i_o); % m

    d_lambda = deg2rad(delta_u + delta_RAAN * cosd(i_o)) * a_o_m; % m
end