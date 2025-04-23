%% AA279D Simulations for PSET 3
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

%%%%% INITIAL RELATIVE ORBITAL ELEMENTS %%%%%%%
d_a_SV2_init = 0; % m
d_lambda_SV2_init = -124350; % m
d_e_x_SV2_init = 110; % m
d_e_y_SV2_init = 202; % m
d_i_x_SV2_init = 79; % m
d_i_y_SV2_init = 1005; % m

d_a_SV3_init = 0; % m
d_lambda_SV3_init = -79328; % m
d_e_x_SV3_init = 42; % m
d_e_y_SV3_init = 452; % m
d_i_x_SV3_init = 36; % m
d_i_y_SV3_init = 827; % m

[r_SV2_init, v_SV2_init] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
d_a_SV2_init,d_lambda_SV2_init,d_e_x_SV2_init,d_e_y_SV2_init,d_i_x_SV2_init,d_i_y_SV2_init);
state_abs_SV2_init = [r_SV2_init;v_SV2_init];

[r_SV3_init, v_SV3_init] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
d_a_SV3_init,d_lambda_SV3_init,d_e_x_SV3_init,d_e_y_SV3_init,d_i_x_SV3_init,d_i_y_SV3_init);
state_abs_SV3_init = [r_SV3_init;v_SV3_init];

% Converting initial swarm ECI coordinates to initial relative co-ordinates in
% chief's RTN frame.
[r_SV2_RTN_init, v_SV2_RTN_init] = ECI2RTN_rel(r_SV1_ECI_init', v_SV1_ECI_init', r_SV2_init', v_SV2_init');
r_SV2_RTN_init = r_SV2_RTN_init';
v_SV2_RTN_init = v_SV2_RTN_init';

[r_SV3_RTN_init, v_SV3_RTN_init] = ECI2RTN_rel(r_SV1_ECI_init', v_SV1_ECI_init', r_SV3_init', v_SV3_init');
r_SV3_RTN_init = r_SV3_RTN_init';
v_SV3_RTN_init = v_SV3_RTN_init';

% This needs to also have the chief's state to propagate
state_rel_SV2_init = [r_SV1_ECI_init; v_SV1_ECI_init; r_SV2_RTN_init; v_SV2_RTN_init];
state_rel_SV3_init = [r_SV1_ECI_init; v_SV1_ECI_init; r_SV3_RTN_init; v_SV3_RTN_init];

%%%%%% RUN SIM OF RELATIVE MOTION %%%%%%%%
[t_3, state3] = rk4_eom_rel_RTN(tstart:tint:tend, state_rel_SV2_init);
[t_4, state4] = rk4_eom_rel_RTN(tstart:tint:tend, state_rel_SV3_init);
SV2_rel_pos = [state3(:, 7), state3(:, 8), state3(:, 9)];
SV2_rel_vel = [state3(:, 10), state3(:, 11), state3(:, 12)];
SV3_rel_pos = [state4(:, 7), state4(:, 8), state4(:, 9)];
SV3_rel_vel = [state4(:, 10), state4(:, 11), state4(:, 12)];

%%%%% RUN ABSOLUTE POSITION SIM WITH SV2 and SV3 %%%%%%
[t_5, state5] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV2_init, false);
[t_6, state6] = rk4_eom_ECI(tstart:tint:tend, state_abs_SV3_init, false);
r_SV2_ECI_no_j2 = state5(:,1:3);
v_SV2_ECI_no_j2 = state5(:,4:6);
r_SV3_ECI_no_j2 = state6(:,1:3);
v_SV3_ECI_no_j2 = state6(:,4:6);

% Convert the absolute ECI positions of the chief and deputies and convert to
% RTN.
[rho_SV2_RTN, rho_SV2_RTN_dot] = ECI2RTN_rel(r_ECI_no_j2, v_ECI_no_j2, r_SV2_ECI_no_j2, v_SV2_ECI_no_j2);
[rho_SV3_RTN, rho_SV3_RTN_dot] = ECI2RTN_rel(r_ECI_no_j2, v_ECI_no_j2, r_SV3_ECI_no_j2, v_SV3_ECI_no_j2);


%%%%% NEW INITIAL RELATIVE ORBIT ELEMENTS %%%%%

d_a_SV2_init_new = 200; % m

d_a_SV3_init_new = -1000; % m

[r_SV2_init, v_SV2_init] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
d_a_SV2_init_new,d_lambda_SV2_init,d_e_x_SV2_init,d_e_y_SV2_init,d_i_x_SV2_init,d_i_y_SV2_init);
state_abs_SV2_init_new = [r_SV2_init;v_SV2_init];

[r_SV3_init, v_SV3_init] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
d_a_SV3_init_new,d_lambda_SV3_init,d_e_x_SV3_init,d_e_y_SV3_init,d_i_x_SV3_init,d_i_y_SV3_init);
state_abs_SV3_init_new = [r_SV3_init;v_SV3_init];


%%%%% APPLY SEMI-MAJOR AXIS MANEUVERs %%%%%
maneuver_time = tend/2; % set arbitrary desired maneuevr time 
% (later find periapsis for fuel-optimality)

desired_delta_a_SV2 = -d_a_SV2_init_new/1000; % m --> km
desired_delta_a_SV3 = -d_a_SV3_init_new/1000; % m --> km

[rho_SV2_RTN_w_maneuver, rho_SV2_RTN_dot_w_maneuver, t_SV2_combined] = ...
apply_sma_maneuver(state_abs_SV2_init_new,tstart,tint,maneuver_time,tend,desired_delta_a_SV2,r_ECI_no_j2,v_ECI_no_j2);

[rho_SV3_RTN_w_maneuver, rho_SV3_RTN_dot_w_maneuver, t_SV3_combined] = ...
apply_sma_maneuver(state_abs_SV3_init_new,tstart,tint,maneuver_time,tend,desired_delta_a_SV3,r_ECI_no_j2,v_ECI_no_j2);


