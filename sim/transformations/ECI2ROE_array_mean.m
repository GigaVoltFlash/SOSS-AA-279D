function [d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y] = ECI2ROE_array_mean(r_o, v_o, r_t, v_t, J2_flag)
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
        % Chief orbit elements, osculating --> mean
        [a_o, e_o, i_o, RAAN_o, omega_o, nu_o, M_o] = ECI2OE(r_o(k,:), v_o(k,:));
        
        osc_elem_o = [a_o*1e3; e_o; deg2rad(i_o); deg2rad(RAAN_o); deg2rad(omega_o); deg2rad(M_o)];
        mean_elem_o = osc2mean(osc_elem_o, J2_flag); 

        a_o = mean_elem_o(1)/1e3; % m --> km
        e_o = mean_elem_o(2);
        i_o = rad2deg(mean_elem_o(3));
        RAAN_o = rad2deg(mean_elem_o(4));
        omega_o = rad2deg(mean_elem_o(5));
        M_o = rad2deg(mean_elem_o(6));

        [~, e_x_o, e_y_o, i_o, RAAN_o, u_o] = OE2quasi_nonsing(a_o, e_o, i_o, RAAN_o, omega_o, M_o);

        % Deputy (target) orbit elements, oscualting --> mean
        [a_t, e_t, i_t, RAAN_t, omega_t, nu_t, M_t] = ECI2OE(r_t(k,:), v_t(k,:));

        osc_elem_t = [a_t*1e3; e_t; deg2rad(i_t); deg2rad(RAAN_t); deg2rad(omega_t); deg2rad(M_t)];
        mean_elem_t = osc2mean(osc_elem_t, J2_flag); 
        
        % Step 6: Convert back to degrees
        a_t = mean_elem_t(1)/1e3; % m --> km
        e_t = mean_elem_t(2);
        i_t = rad2deg(mean_elem_t(3));
        RAAN_t = rad2deg(mean_elem_t(4));
        omega_t = rad2deg(mean_elem_t(5));
        M_t = rad2deg(mean_elem_t(6));

        [~, e_x_t, e_y_t, i_t, RAAN_t, u_t] = OE2quasi_nonsing(a_t, e_t, i_t, RAAN_t, omega_t, M_t);

        % Relative orbital elements (ROE)
        [d_a(k), d_lambda(k), d_e_x(k), d_e_y(k), d_i_x(k), d_i_y(k)] = ...
            quasi_nonsing2ROE(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o, ...
                              a_t, e_x_t, e_y_t, i_t, RAAN_t, u_t);
    end
end