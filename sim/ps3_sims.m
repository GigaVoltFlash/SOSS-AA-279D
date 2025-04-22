%% AA279D Simulations for PSET 3
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

% Propagate the chief to get r0 for future analysis
ps1_sims;

% Define initial relative orbital elements for Q1 (small separation,
% minimal eccentricity)
% I reduced the lambda value by a few orders of magnitude, and took the
% eccentricity to zero for both spacecraft

d_a_SV2_init = 0; % m
d_lambda_SV2_init = 0; % m % -124
d_e_x_SV2_init = 0; % m
d_e_y_SV2_init = 100; % m
d_i_x_SV2_init = 0; % m % 79
d_i_y_SV2_init = 1000; % m % 1005

d_a_SV3_init = 0; % m
d_lambda_SV3_init = 0; % m % -793
d_e_x_SV3_init = 0; % m
d_e_y_SV3_init = 200; % m
d_i_x_SV3_init = 0; % m
d_i_y_SV3_init = 800; % m % 827

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
[t_2, state_2] = rk4_eom_rel_RTN(tstart:tint:tend, state_rel_SV2_init);
[t_3, state_3] = rk4_eom_rel_RTN(tstart:tint:tend, state_rel_SV3_init);
SV2_rel_pos = [state_2(:, 7), state_2(:, 8), state_2(:, 9)];
SV2_rel_vel = [state_2(:, 10), state_2(:, 11), state_2(:, 12)];
SV3_rel_pos = [state_3(:, 7), state_3(:, 8), state_3(:, 9)];
SV3_rel_vel = [state_3(:, 10), state_3(:, 11), state_3(:, 12)];

% We see that the ratio of relative distance to orbit size is smaller than 1e-3
rho_pos_ratio_SV2 = vecnorm(SV2_rel_pos, 2, 2)./r_RTN_no_j2(:, 1);
rho_pos_ratio_SV3 = vecnorm(SV3_rel_pos, 2, 2)./r_RTN_no_j2(:, 1);

%%%%% CALCULATE THE HCW INTEGRATION CONSTANTS %%%%% 
K_SV2 = hcw_integration_constants([r_SV2_RTN_init; v_SV2_RTN_init], a_SV1_init);
K_SV3 = hcw_integration_constants([r_SV3_RTN_init; v_SV3_RTN_init], a_SV1_init);

%%%%% EVALUATE HCW (TWO DIFFERENT METHODS %%%%%%%%%
SV2_HCW_state = evaluate_HCW(t_2, a_SV1_init, K_SV2);
SV3_HCW_state = evaluate_HCW(t_2, a_SV1_init, K_SV3);
SV2_HCW_pos_test = hcw_equations_with_ic(t_2, a_SV1_init, [r_SV2_RTN_init; v_SV2_RTN_init]);
SV3_HCW_pos_test = hcw_equations_with_ic(t_2, a_SV1_init, [r_SV3_RTN_init; v_SV3_RTN_init]);

SV2_HCW_pos = SV2_HCW_state(:, 1:3);
SV3_HCW_pos = SV3_HCW_state(:, 1:3);
SV2_HCW_vel = SV2_HCW_state(:, 4:6);
SV3_HCW_vel = SV3_HCW_state(:, 4:6);



%%%%%%% QUESTION 2 %%%%%%%%%% Need to clean up variables because just
%%%%%%% repeating Q1 rn

%%%%%%% DEFINE NEW ECCENTRIC INITIAL CONDITIONS %%%%%%%
% Define initial relative orbital elements for Q2 (small separation,
% significant eccentricty)
% Same as Q1, except eccentricity is now 0.1

load_SV1_eccentric;


