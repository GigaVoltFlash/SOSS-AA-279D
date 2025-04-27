function [d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y] = ECI2ROE_array(r_o, v_o, r_t, v_t)
% Takes in arrays of position and velocity of the chief and deputy, 
% converts to quasi-nonsingular ROE of deputy

    N = size(r_o, 1);

    d_a = zeros(N,1);
    d_lambda = zeros(N,1);
    d_e_x = zeros(N,1);
    d_e_y = zeros(N,1);
    d_i_x = zeros(N,1);
    d_i_y = zeros(N,1);

    for k = 1:N
        % Chief orbit elements
        [a_o, e_o, i_o, RAAN_o, omega_o, nu_o, M_o] = ECI2OE(r_o(k,:), v_o(k,:));
        [~, e_x_o, e_y_o, i_o, RAAN_o, u_o] = OE2quasi_nonsing(a_o, e_o, i_o, RAAN_o, omega_o, M_o);


        % Deputy (target) orbit elements
        [a_t, e_t, i_t, RAAN_t, omega_t, nu_t, M_t] = ECI2OE(r_t(k,:), v_t(k,:));
        [~, e_x_t, e_y_t, i_t, RAAN_t, u_t] = OE2quasi_nonsing(a_t, e_t, i_t, RAAN_t, omega_t, M_t);

        % Relative orbital elements (ROE)
        [d_a(k), d_lambda(k), d_e_x(k), d_e_y(k), d_i_x(k), d_i_y(k)] = ...
            quasi_nonsing2ROE(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o, ...
                              a_t, e_x_t, e_y_t, i_t, RAAN_t, u_t);
    end
end