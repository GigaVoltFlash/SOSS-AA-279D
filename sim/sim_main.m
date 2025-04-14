%% AA279D Master Simulation
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 
close all;
addpath("plotting_methods\");
addpath("mean_osc\");

% Load all constants from the constants file
constants()
global mu_earth

% Flags for which plots to create
PS1_plots = false;
PS2_plots = true;

% Units: kilometers and degrees
a_SV1_init = 6944;
ex_SV1_init = -4e-5;
ey_SV1_init = 1.6e-3;
i_SV1_init = 99.4;
RAAN_SV1_init = -151.1;
u_SV1_init = -47.9;
[a_SV1_init, e_SV1_init, i_SV1_init, RAAN_SV1_init, w_SV1_init, nu_SV1_init, M_SV1_init] = quasi_nonsing2OE(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init);

[r_SV1_ECI_init, v_SV1_ECI_init] = OE2ECI(a_SV1_init, e_SV1_init, i_SV1_init, RAAN_SV1_init, w_SV1_init, nu_SV1_init);
state_init = [r_SV1_ECI_init; v_SV1_ECI_init];

%%% PROPAGATION TIME DETAILS %%%
tstart = 0; % Start time (s)
t_orbit = 2*pi*sqrt(a_SV1_init^3/mu_earth);
tend = 25*t_orbit; % End time (s)
tint = t_orbit/500.0; % Time step (s)

%%%%%%%%%%% RUN SIM WITH AND WITHOUT J2 %%%%%%%%%%%
[t_1, state1] = rk4_eom_ECI(tstart:tint:tend, state_init, false);
[t_2, state2] = rk4_eom_ECI(tstart:tint:tend, state_init, true);
r_ECI_with_j2 = state2(:,1:3);
v_ECI_with_j2 = state2(:,4:6);
r_ECI_no_j2 = state1(:,1:3);
v_ECI_no_j2 = state1(:,4:6);

%%%%%%% KEPLERIAN PROPAGATION %%%%%%% 
% Propagates over the same times t_1 as the sim without J2 
[r_ECI_keplerian, v_ECI_keplerian] = keplerian_propagator(a_SV1_init, e_SV1_init, i_SV1_init, RAAN_SV1_init, w_SV1_init, t_1, M_SV1_init);

%%%%%% CONVERT TO RTN %%%%%%%%
[r_RTN_no_j2, v_RTN_no_j2] = ECI2RTN(r_ECI_no_j2,v_ECI_no_j2);
[r_RTN_keplerian, v_RTN_keplerian] = ECI2RTN(r_ECI_keplerian,v_ECI_keplerian);

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
% % % d_a_2_init_new = 1000; % m
% % % 
% % % % redo b and c: states 7,8,9,10
% % % 
% % % 
% % % 
% % % %%%%% APPLY SEMI-MAJOR AXIS MANEUVER %%%%%
% % % desired_delta_a = 12; % km
% % % maneuver_time = tend/2;
% % % 
% % % [t_11, state11] = rk4_eom_ECI(tstart:tint:maneuver_time, state_abs_SV2_init, false);
% % % r_SV2_ECI_no_j2_inter = state11(end,1:3);
% % % v_SV2_ECI_no_j2_inter = state11(end,4:6);
% % % [a_SV2_inter,e_SV2_inter] = ECI2OE(r_SV2_ECI_no_j2_inter,v_SV2_ECI_no_j2_inter);
% % % 
% % % delta_v_t_RTN = sma_maneuver(a_SV2_inter,e_SV2_inter,r_SV2_ECI_no_j2_inter,desired_delta_a);
% % % 
% % % delta_v_ECI = dv_RTN2ECI(r_SV2_ECI_no_j2_inter, v_SV2_ECI_no_j2_inter, delta_v_t_RTN);
% % % % do for SV3, state 12
% % % v_SV2_ECI_no_j2_new = v_SV2_ECI_no_j2_inter+delta_v_ECI;
% % % state_abs_SV2_inter_new = [r_SV2_ECI_no_j2_inter;v_SV2_ECI_no_j2_inter_new];
% % % 
% % % [t_13, state13] = rk4_eom_ECI(maneuver_time:tint:tend, state_abs_SV2_inter_new, false);
% % % % do for SV3, state 12
% % % 
% % % % Plot result showing bounded relative motion


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if PS1_plots
    %%%%%% PLOT 3D COMPARISON PLOTS %%%%%% 
    plot_3D_orbit_aoe_compare(r_ECI_no_j2, r_ECI_with_j2, 'Orbit without J2', 'Orbit with J2')
    plot_3D_orbit_aoe_compare(r_ECI_no_j2, r_ECI_keplerian, 'Orbit without J2', 'Keplerian Propagator')
    
    %%%%%% COMPARE RTN VALUES %%%%%% 
    plot_rtn_compare(t_1, r_RTN_no_j2, r_RTN_keplerian, v_RTN_no_j2, v_RTN_keplerian, 'Numerical with no J2', 'Keplerian', 'Error between No-J2 Numerical and Keplerian, large time-step');
    
    %%%%%%% COMPUTE AND PLOT OE, ECCENTRICITY, ANGULAR MOMENTUM, SPECIFIC MECHANICAL ENERGY %%%%%%% 
    [a_2,e_2,i_2,RAAN_2,omega_2,nu_2] = compute_and_plot_orbit_params(r_ECI_no_j2,v_ECI_no_j2,r_ECI_with_j2,v_ECI_with_j2,t_1,t_2);
    
    %%%%%%%% COMPUTE AND PLOT OSCULATING VERSUS MEAN %%%%%%%%%%%%
    [a_3,e_3,i_3,RAAN_3,omega_3,nu_3,t_3] = compute_mean_oe(a_SV1_init,ex_init,ey_init,i_init,RAAN_init,u_init,tstart,tint,tend);
    plot_osc_mean_oe(a_2,e_2,i_2,RAAN_2,omega_2,nu_2,t_2,a_3,e_3,i_3,RAAN_3,omega_3,nu_3,t_3);
end

if PS2_plots

    %%%%%% PLOT OF ABSOLUTE ORBITS OF SV1 AND SV2 %%%%%%%%%%
    plot_3D_orbits(r_ECI_no_j2, r_SV2_ECI_no_j2, 'Chief', 'Deputy')
    
    %%%%%%%% PLOT RELATIVE POSITIONS (IN CHIEF'S RTN FRAME) %%%%%%%
    plot_rel_sat_pos_3D(SV2_rel_pos, 'figures/SV2_3d_traj_rel_eom.png') % For part b
    plot_rel_sat_pos_3D(rho_SV2_RTN, 'figures/SV2_3d_traj_eci2rtn.png') % For part c
    plot_rel_pos_vel(t_5, rho_SV2_RTN, rho_SV2_RTN_dot, 'SV2 Relative RTN position and velocity using ECI interpolation', 'figures/SV2_rel_pos_vel_eci2rtn.png'); % For part c
    plot_rel_pos_vel(t_5, rho_SV2_RTN, rho_SV2_RTN_dot, 'SV2 Relative RTN position and velocity using non-linear relative EOMs', 'figures/SV2_rel_pos_vel_eom.png'); % For part b

    plot_rel_sat_pos_3D(SV3_rel_pos, 'figures/SV3_3d_traj_rel_eom.png') % For part b
    plot_rel_sat_pos_3D(rho_SV3_RTN, 'figures/SV3_3d_traj_eci2rtn.png') % For part c
    plot_rel_pos_vel(t_6, rho_SV3_RTN, rho_SV3_RTN_dot, 'SV3 Relative RTN position and velocity using ECI interpolation', 'figures/SV3_rel_pos_vel_eci2rtn.png'); % For part c
    plot_rel_pos_vel(t_6, rho_SV3_RTN, rho_SV3_RTN_dot, 'SV3 Relative RTN position and velocity using non-linear relative EOMs', 'figures/SV3_rel_pos_vel_eom.png'); % For part b

    plot_RT_RN_projections(SV2_rel_pos, SV3_rel_pos, 'figures/RTN_projection_rel_eom.png') % For part b
    plot_RT_RN_projections(rho_SV2_RTN, rho_SV3_RTN, 'figures/RTN_projection_eci2rtn.png') % For part c

    compare_rel_pos_error(t_2, rho_SV2_RTN, SV2_rel_pos, rho_SV2_RTN_dot, SV2_rel_vel, 'figures/SV2_error_in_rel_methods.png'); % For part d
    compare_rel_pos_error(t_2, rho_SV3_RTN, SV3_rel_pos, rho_SV3_RTN_dot, SV3_rel_vel, 'figures/SV3_error_in_rel_methods.png'); % For part d

end

