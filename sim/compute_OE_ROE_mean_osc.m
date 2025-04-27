function [d_a_osc, d_lambda_osc, d_e_x_osc, d_e_y_osc, d_i_x_osc, d_i_y_osc, ...
          d_a_mean, d_lambda_mean, d_e_x_mean, d_e_y_mean, d_i_x_mean, d_i_y_mean, ...
          a_osc, e_x_osc, e_y_osc, i_osc, RAAN_osc, u_osc, ...
          a_mean, e_x_mean, e_y_mean, i_mean, RAAN_mean, u_mean] = compute_OE_ROE_mean_osc(r_o, v_o, r_t, v_t, J2_flag)
    % Inputs:
    %   r_o, v_o : chief position and velocity (N x 3)
    %   r_t, v_t : deputy position and velocity (N x 3)
    %   J2_flag  : 1 = apply J2 corrections, 0 = no J2
    %
    % Outputs:
    %   Osculating ROEs: (d_a_osc, d_lambda_osc, d_e_x_osc, d_e_y_osc, d_i_x_osc, d_i_y_osc)
    %   Mean ROEs: (d_a_mean, d_lambda_mean, d_e_x_mean, d_e_y_mean, d_i_x_mean, d_i_y_mean)
    %   Osculating quasi-nonsingular elements (a_osc, e_x_osc, e_y_osc, i_osc, RAAN_osc, u_osc)
    %   Mean quasi-nonsingular elements (a_mean, e_x_mean, e_y_mean, i_mean, RAAN_mean, u_mean)

    % Osculating quasi-nonsingular orbital elements
    [a_osc,e_osc,i_osc,RAAN_osc,omega_osc,nu_osc,M_osc] = ECI2OE_array(r_t,v_t); % degree outputs
    [a_osc,e_x_osc,e_y_osc,i_osc,RAAN_osc,u_osc] = OE2quasi_nonsing_array(a_osc,e_osc,i_osc,RAAN_osc,omega_osc,M_osc); % degree inputs and outputs
    
    % Mean quasi-non-singular orbital elements    
    osc_elem_array = [a_osc*1e3; e_osc; deg2rad(i_osc); deg2rad(RAAN_osc); deg2rad(omega_osc); deg2rad(M_osc)]; % 6 x N
    
    mean_elem_array = zeros(size(osc_elem_array)); % Preallocate
    for k = 1:size(osc_elem_array, 2)
        mean_elem_array(:,k) = osc2mean(osc_elem_array(:,k), J2_flag);
    end
    
    a_mean = mean_elem_array(1,:) / 1e3; % m --> km
    e_mean = mean_elem_array(2,:);
    i_mean_deg = rad2deg(mean_elem_array(3,:));
    RAAN_mean_deg = rad2deg(mean_elem_array(4,:));
    omega_mean_deg = rad2deg(mean_elem_array(5,:));
    M_mean_deg = rad2deg(mean_elem_array(6,:));
    
    [a_mean, e_x_mean, e_y_mean, i_mean, RAAN_mean, u_mean] = ...
        OE2quasi_nonsing_array(a_mean, e_mean, i_mean_deg, RAAN_mean_deg, omega_mean_deg, M_mean_deg); % deg inputs and outputs
    
    % Osculating quasi-nonsingular relative orbital elements
    [d_a_osc, d_lambda_osc, d_e_x_osc, d_e_y_osc, d_i_x_osc, d_i_y_osc] = ...
        ECI2ROE_array(r_o, v_o, r_t, v_t);
    
    % Mean quasi-nonsingular relative orbital elements
    [d_a_mean, d_lambda_mean, d_e_x_mean, d_e_y_mean, d_i_x_mean, d_i_y_mean] = ...
        ECI2ROE_array_mean(r_o, v_o, r_t, v_t, J2_flag);
end