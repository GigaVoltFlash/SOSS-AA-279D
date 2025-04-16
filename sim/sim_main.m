%% AA279D Master Simulation
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 
close all;
addpath("plotting_methods\");
addpath("mean_osc\");
addpath("transformations\");
addpath("propagators\");
addpath("eoms\");

% Load all constants from the constants file
constants()
global mu_earth;

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

d_a_SV2_init_new = 200; % m

d_a_SV3_init_new = -1000; % m

[r_SV2_init, v_SV2_init] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
d_a_SV2_init_new,d_lambda_SV2_init,d_e_x_SV2_init,d_e_y_SV2_init,d_i_x_SV2_init,d_i_y_SV2_init);
state_abs_SV2_init_new = [r_SV2_init;v_SV2_init];

[r_SV3_init, v_SV3_init] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
d_a_SV3_init_new,d_lambda_SV3_init,d_e_x_SV3_init,d_e_y_SV3_init,d_i_x_SV3_init,d_i_y_SV3_init);
state_abs_SV3_init_new = [r_SV3_init;v_SV3_init];


% %%%%% APPLY SEMI-MAJOR AXIS MANEUVERs %%%%%
% maneuver_time = tend/2; % set arbitrary desired maneuevr time 
% % (later find periapsis for fuel-optimality)
% 
% desired_delta_a_SV2 = -d_a_SV2_init_new/1000; % m --> km
% desired_delta_a_SV3 = -d_a_SV3_init_new/1000; % m --> km
% 
% [rho_SV2_RTN_w_maneuver, rho_SV2_RTN_dot_w_maneuver, t_SV2_combined] = ...
% apply_sma_maneuver(state_abs_SV2_init_new,tstart,tint,maneuver_time,tend,desired_delta_a_SV2,r_ECI_no_j2,v_ECI_no_j2);
% 
% [rho_SV3_RTN_w_maneuver, rho_SV3_RTN_dot_w_maneuver, t_SV3_combined] = ...
% apply_sma_maneuver(state_abs_SV3_init_new,tstart,tint,maneuver_time,tend,desired_delta_a_SV3,r_ECI_no_j2,v_ECI_no_j2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if PS1_plots
    %%%%%% PLOT 3D COMPARISON PLOTS %%%%%% 
    plot_3D_orbit_aoe_compare(r_ECI_no_j2, r_ECI_with_j2, 'Orbit without J2', 'Orbit with J2')
    plot_3D_orbit_aoe_compare(r_ECI_no_j2, r_ECI_keplerian, 'Orbit without J2', 'Keplerian Propagator')
    
    %%%%%% COMPARE RTN VALUES %%%%%% 
    plot_rtn_compare(t_1, t_orbit, r_RTN_no_j2, r_RTN_keplerian, v_RTN_no_j2, v_RTN_keplerian, 'Numerical with no J2', 'Keplerian', 'Error between No-J2 Numerical and Keplerian, small time-step', 'figures/comparing_rtn_small_timestep.png');
    
    %%%%%%% COMPUTE AND PLOT OE, ECCENTRICITY, ANGULAR MOMENTUM, SPECIFIC MECHANICAL ENERGY %%%%%%% 
    [a_2,e_2,i_2,RAAN_2,omega_2,nu_2] = compute_and_plot_orbit_params(r_ECI_no_j2,v_ECI_no_j2,r_ECI_with_j2,v_ECI_with_j2,t_1,t_2);
    
    %%%%%%%% COMPUTE AND PLOT OSCULATING VERSUS MEAN %%%%%%%%%%%%
    [a_3,e_3,i_3,RAAN_3,omega_3,nu_3,t_3] = compute_mean_oe(a_SV1_init,ex_init,ey_init,i_init,RAAN_init,u_init,tstart,tint,tend);
    plot_osc_mean_oe(a_2,e_2,i_2,RAAN_2,omega_2,nu_2,t_2,a_3,e_3,i_3,RAAN_3,omega_3,nu_3,t_3);
end

if PS2_plots

    %%%%%% PLOT OF ABSOLUTE ORBITS OF SV1 AND SV2 %%%%%%%%%%
    %plot_3D_orbits(r_ECI_no_j2, r_SV2_ECI_no_j2, 'Chief', 'Deputy')
    
    %%%%%%%% PLOT RELATIVE POSITIONS (IN CHIEF'S RTN FRAME) %%%%%%%
    
    plot_rel_pos_vel(t_3, t_orbit, SV2_rel_pos, SV2_rel_vel, rho_SV2_RTN, rho_SV2_RTN_dot, 'SV2 RTN position and velocity', 'figures/SV2_rel_pos_vel.png'); % For part b and c
    plot_rel_pos_vel(t_4, t_orbit, SV3_rel_pos, SV3_rel_vel, rho_SV3_RTN, rho_SV3_RTN_dot, 'SV3 RTN position and velocity', 'figures/SV3_rel_pos_vel.png'); % For part b and c

    % plot_rel_pos_vel(t_SV2_combined, t_orbit, rho_SV2_RTN_w_maneuver, rho_SV2_RTN_dot_w_maneuver, 'SV2 RTN position and velocity with maneuver','figures/SV2_rel_pos_vel_eom_maneuver.png'); % For part f
    % plot_rel_pos_vel(t_SV3_combined, t_orbit, rho_SV3_RTN_w_maneuver, rho_SV3_RTN_dot_w_maneuver, 'SV3 RTN position and velocity with maneuver','figures/SV3_rel_pos_vel_eom_maneuver.png'); % For part f
    
    %%%%%%% 3D PLOTTING OF ORBITS %%%%%%%%%
    plot_rel_sat_pos_3D_multi_method(SV2_rel_pos, rho_SV2_RTN, SV3_rel_pos, rho_SV3_RTN, 'figures/SV2_SV3_3d_traj_rel.png') % For part b and c
    % plot_rel_sat_pos_3D(rho_SV2_RTN_w_maneuver, 'figures/SV2_3d_traj_rel_eom.png') % For part f
    
    %%%%%%% PROJECTING ONTO RTN PLANES %%%%%%%%%
    plot_RT_RN_projections_both(SV2_rel_pos, rho_SV2_RTN, SV3_rel_pos, rho_SV3_RTN, 'RTN position of SV2 and SV3','figures/RTN_projection.png') % For part c
    
    %%%%%%% COMPARING ERRORS BETWEEN METHODS %%%%%%%%%
    compare_rel_pos_error(t_5, t_orbit, rho_SV2_RTN, SV2_rel_pos, rho_SV2_RTN_dot, SV2_rel_vel,'SV2 RTN Error Comparison','figures/SV2_error_in_rel_methods.png'); % For part d
    compare_rel_pos_error(t_6, t_orbit,  rho_SV3_RTN, SV3_rel_pos, rho_SV3_RTN_dot, SV3_rel_vel,'SV3 RTN Error Comparison','figures/SV3_error_in_rel_methods.png'); % For part d
    % compare_rel_pos_error(t_5, t_orbit, rho_SV2_RTN, SV2_rel_pos, rho_SV2_RTN_dot, SV2_rel_vel,'SV2 RTN Error Comparison with non-zero relative SMA','figures/SV2_error_in_rel_methods_nonzero_a.png'); % For part d, run with non-zero delta a
end

