%% AA279D Simulations for PSET 4
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

%%%%% INITIAL RELATIVE ORBITAL ELEMENTS %%%%%%% 
% From PS3
d_a_SV2_init_1 = 0; % m
d_lambda_SV2_init_1 = 0; % m % -124
d_e_x_SV2_init_1 = 0; % m
d_e_y_SV2_init_1 = 100; % m
d_i_x_SV2_init_1 = 0; % m % 79
d_i_y_SV2_init_1 = 1000; % m % 1005

d_a_SV3_init_1 = 0; % m
d_lambda_SV3_init_1 = 0; % m % -793
d_e_x_SV3_init_1 = 0; % m
d_e_y_SV3_init_1 = 200; % m
d_i_x_SV3_init_1 = 0; % m
d_i_y_SV3_init_1 = 800; % m % 827

% Given in PS4
d_a_SV2_init_2 = 0; % m
d_lambda_SV2_init_2 = 100; % m % -124
d_e_x_SV2_init_2 = 50; % m
d_e_y_SV2_init_2 = 100; % m
d_i_x_SV2_init_2 = 30; % m % 79
d_i_y_SV2_init_2 = 200; % m % 1005

% Applying similar elements for SV3
d_a_SV3_init_2 = 0; % m
d_lambda_SV3_init_2 = -100; % m % -793
d_e_x_SV3_init_2 = 50; % m
d_e_y_SV3_init_2 = 100; % m
d_i_x_SV3_init_2 = 30; % m
d_i_y_SV3_init_2 = 200; % m % 827

% Converting to ECI position and velocities (initial state)
[r_SV2_init_1, v_SV2_init_1] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
d_a_SV2_init_1,d_lambda_SV2_init_1,d_e_x_SV2_init_1,d_e_y_SV2_init_1,d_i_x_SV2_init_1,d_i_y_SV2_init_1);
state_abs_SV2_init_1 = [r_SV2_init_1;v_SV2_init_1];

[r_SV3_init_1, v_SV3_init_1] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
d_a_SV3_init_1,d_lambda_SV3_init_1,d_e_x_SV3_init_1,d_e_y_SV3_init_1,d_i_x_SV3_init_1,d_i_y_SV3_init_1);
state_abs_SV3_init_1 = [r_SV3_init_1;v_SV3_init_1];

[r_SV2_init_2, v_SV2_init_2] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
d_a_SV2_init_2,d_lambda_SV2_init_2,d_e_x_SV2_init_2,d_e_y_SV2_init_2,d_i_x_SV2_init_2,d_i_y_SV2_init_2);
state_abs_SV2_init_2 = [r_SV2_init_2;v_SV2_init_2];

[r_SV3_init_2, v_SV3_init_2] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
d_a_SV3_init_2,d_lambda_SV3_init_2,d_e_x_SV3_init_2,d_e_y_SV3_init_2,d_i_x_SV3_init_2,d_i_y_SV3_init_2);
state_abs_SV3_init_2 = [r_SV3_init_2;v_SV3_init_2];

% Chief
state_abs_SV1_init = [r_SV1_ECI_init, v_SV1_ECI_init];

%%%%% RUN ABSOLUTE POSITION SIM WITH SV2 and SV3 %%%%%%
% Chief no J2
[t_3, state3] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV1_init, false);
r_SV1_no_j2 = state3(:,1:3);
v_SV1_no_j2 = state3(:,4:6);

% Chief with J2
[t_4, state4] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV1_init, true);
r_SV1_with_j2 = state4(:,1:3);
v_SV1_with_j2 = state4(:,4:6);

% Initial conditions 1 and no J2
[t_5, state5] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV2_init_1, false);
[t_6, state6] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV3_init_1, false);
r_SV2_1_no_j2 = state5(:,1:3);
v_SV2_1_no_j2 = state5(:,4:6);
r_SV3_1_no_j2 = state6(:,1:3);
v_SV3_1_no_j2 = state6(:,4:6);

% Initial conditions 1 and with J2
[t_7, state7] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV2_init_1, true);
[t_8, state8] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV3_init_1, true);
r_SV2_1_with_j2 = state7(:,1:3);
v_SV2_1_with_j2 = state7(:,4:6);
r_SV3_1_with_j2 = state8(:,1:3);
v_SV3_1_with_j2 = state8(:,4:6);

% Initial conditions 2 and no J2
[t_9, state9] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV2_init_2, false);
[t_10, state10] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV3_init_2, false);
r_SV2_2_no_j2 = state9(:,1:3);
v_SV2_2_no_j2 = state9(:,4:6);
r_SV3_2_no_j2 = state10(:,1:3);
v_SV3_2_no_j2 = state10(:,4:6);

% Initial conditions 1 and with J2
[t_11, state11] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV2_init_2, true);
[t_12, state12] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV3_init_2, true);
r_SV2_2_with_j2 = state11(:,1:3);
v_SV2_2_with_j2 = state11(:,4:6);
r_SV3_2_with_j2 = state12(:,1:3);
v_SV3_2_with_j2 = state12(:,4:6);

%%%%%%% CALCULATE OSC AND MEAN OE %%%%%%%%%%%

% Osculating quasi-nonsingular orbital elements
[a_osc,e_osc,i_osc,RAAN_osc,omega_osc,nu_osc,M_osc] = ECI2OE_array(r,v); % degree outputs
[a_osc,e_x_osc,e_y_osc,i_osc,RAAN_osc,u_osc] = OE2quasi_nonsing_array(a_osc,e_osc,i_osc,RAAN_osc,omega_osc,M_osc); % degree inputs and outputs

% Mean quasi-non-singular orbital elements
J2_flag = 0;% 1 = J2, 0 = no J2
i_osc_rad = deg2rad(i_osc); RAAN_osc_rad = deg2rad(RAAN_osc);
omega_osc_rad = deg2rad(omega_osc); M_osc_rad = deg2rad(M_osc);
[a_mean,e_mean,i_mean_rad,RAAN_mean_rad,omega_mean_rad,M_mean_rad] = ...
    osc2mean(a_osc,e_osc,i_osc_rad,RAAN_osc_rad,omega_osc_rad,M_osc_rad,J2_flag); % rad inputs and ouputs
i_mean_deg = rad2deg(i_mean_rad); RAAN_mean_deg = rad2deg(RAAN_mean_rad);
omega_mean_deg = rad2deg(omega_mean_rad); M_mean_deg = rad2deg(M_mean_rad);
[a_mean,e_x_mean,e_y_mean,i_mean,RAAN_mean,u_mean] = OE2quasi_nonsing_array(a_mean,e_mean,i_mean_deg,RAAN_mean_deg,omega_mean_deg,M_mean_deg); % degree inputs and outputs





%%%%%%%%%%%%% CONVERT TO RTN %%%%%%%%%%%%%%%%
% Convert the absolute ECI positions of the chief and deputies and convert to RTN.
%[rho_SV2_RTN, rho_SV2_RTN_dot] = ECI2RTN_rel(r_SV1_ECI_no_j2, v_SV1_ECI_no_j2, r_SV2_ECI_no_j2, v_SV2_ECI_no_j2);
%[rho_SV3_RTN, rho_SV3_RTN_dot] = ECI2RTN_rel(r_SV1_ECI_no_j2, v_SV1_ECI_no_j2, r_SV3_ECI_no_j2, v_SV3_ECI_no_j2);