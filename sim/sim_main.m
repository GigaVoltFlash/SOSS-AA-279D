%% AA279D Master Simulation
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 
close all;
clear;
addpath("plotting_methods\");
addpath("mean_osc\");
addpath("transformations\");
addpath("propagators\");
addpath("eoms\");

% Load all constants from the constants file
constants;

% Flags for which plots to create
run_PS1 = false;
run_PS2 = false;
run_PS3 = true;

% Load SV1 information
load_SV1;

%%% PROPAGATION TIME DETAILS %%%
tstart = 0; % Start time (s)
t_orbit = 2*pi*sqrt(a_SV1_init^3/mu_earth);
tend = 25*t_orbit; % End time (s)
tint = t_orbit/500.0; % Time step (s)

% Choose which simulations to run
if run_PS1
    ps1_sims;
end
if run_PS2
    ps2_sims;
end
if run_PS3
    ps3_sims;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if run_PS1
    %%%%%% PLOT 3D COMPARISON PLOTS %%%%%% 
    plot_3D_orbit_aoe_compare(r_ECI_no_j2, r_ECI_with_j2, 'Orbit without J2', 'Orbit with J2')
    plot_3D_orbit_aoe_compare(r_ECI_no_j2, r_ECI_keplerian, 'Orbit without J2', 'Keplerian Propagator')
    
    %%%%%% COMPARE RTN VALUES %%%%%% 
    % plot_rtn_compare(t_1, t_orbit, r_RTN_no_j2, r_RTN_keplerian, v_RTN_no_j2, v_RTN_keplerian, 'Numerical with no J2', 'Keplerian', 'Error between No-J2 Numerical and Keplerian, large time-step', 'figures/comparing_rtn_large_timestep.png');
    plot_rtn_compare(t_1, t_orbit, r_RTN_no_j2, r_RTN_keplerian, v_RTN_no_j2, v_RTN_keplerian, 'Numerical with no J2', 'Keplerian', 'Error between No-J2 Numerical and Keplerian, small time-step', 'figures/comparing_rtn_small_timestep.png');
    
    %%%%%%% COMPUTE AND PLOT OE, ECCENTRICITY, ANGULAR MOMENTUM, SPECIFIC MECHANICAL ENERGY %%%%%%% 
    [a_2,e_2,i_2,RAAN_2,omega_2,nu_2] = compute_and_plot_orbit_params(r_ECI_no_j2,v_ECI_no_j2,r_ECI_with_j2,v_ECI_with_j2,t_1,t_2);
    
    %%%%%%%% COMPUTE AND PLOT OSCULATING VERSUS MEAN %%%%%%%%%%%%
    [a_3,e_3,i_3,RAAN_3,omega_3,nu_3,t_3] = compute_mean_oe(a_SV1_init,ex_init,ey_init,i_init,RAAN_init,u_init,tstart,tint,tend);
    plot_osc_mean_oe(a_2,e_2,i_2,RAAN_2,omega_2,nu_2,t_2,a_3,e_3,i_3,RAAN_3,omega_3,nu_3,t_3);
end

if run_PS2

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

if run_PS3
    %%%%% PLOT OF RELATIVE ORBITS USING HCW %%%%%
    plot_RT_RN_projections(SV2_HCW_pos, SV3_HCW_pos, 'HCW evaluated relative orbits', 'figures/PS3/RTN_projections_HCW.png');
    plot_RT_RN_projections(SV2_rel_pos, SV3_rel_pos, 'HCW evaluated relative orbits', 'figures/PS3/RTN_projections_numerical.png');
end
