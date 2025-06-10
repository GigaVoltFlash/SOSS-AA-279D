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
addpath(genpath("..\gnc\"));

% Load all constants from the constants file
constants;

% Flags for which plots to create
run_PS1 = false;
run_PS2 = false;
run_PS3 = true;
run_PS4 = false;
run_PS5 = false;
run_PS6 = false; 
run_PS7 = false;
run_PS8 = false;
run_PS9 = false;

% Load SV1 information
load_SV1;

%%% PROPAGATION TIME DETAILS %%%
num_orbits = 30; 
tstart = 0; % Start time (s)
t_orbit = 2*pi*sqrt(a_SV1_init^3/mu_earth);
tend = num_orbits*t_orbit; % End time (s)
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
if run_PS4
    ps4_sims;
end
if run_PS5
    ps5_sims;
end
if run_PS6
    ps6_sims;
end
if run_PS7
    ps7_sims;
end
if run_PS8
    ps8_sims;
end
if run_PS9
    ps9_sims;
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
    %plot_rel_pos_vel_single(t_2, t_orbit, SV2_HCW_pos, SV2_HCW_vel, 'SV2 HCW evaluate relative position + velocity', 'figures/PS3/HCW_pos_vel_SV2.png');
    %plot_rel_pos_vel_single(t_2, t_orbit, SV3_HCW_pos, SV3_HCW_vel, 'SV3 HCW evaluate relative position + velocity', 'figures/PS3/HCW_pos_vel_SV3.png');
    %plot_RT_RN_projections(SV2_HCW_pos_test, SV3_HCW_pos_test, 'HCW evaluated relative orbits', 'figures/PS3/RTN_projections_HCW.png');
    %plot_RT_RN_projections(SV2_rel_pos, SV3_rel_pos, 'Non-linear EOMs relative orbits', 'figures/PS3/RTN_projections_numerical.png');

    %plot_RT_RN_projections_both(SV2_HCW_pos,SV2_rel_pos,SV3_HCW_pos,SV3_rel_pos, 'Testing HCW comparison', 'figures/PS3/RTN_projections_HCW_comparison.png');
    %plot_3D_rel_orbit(SV2_HCW_pos,SV3_HCW_pos, 'figures/PS3/3D_HCW_comparison.png');
    %plot_rel_pos_vel_single(t_2, t_orbit, SV2_HCW_pos, SV2_HCW_vel, 'SV2 HCW evaluate relative position + velocity', 'figures/PS3/HCW_pos_vel_SV2.png');
    %plot_rel_pos_vel_single(t_2, t_orbit, SV3_HCW_pos, SV3_HCW_vel, 'SV3 HCW evaluate relative position + velocity', 'figures/PS3/HCW_pos_vel_SV3.png');
    %plot_RT_RN_projections(SV2_HCW_pos, SV3_HCW_pos, 'HCW evaluated relative orbits', 'figures/PS3/RTN_projections_HCW.png');
    % plot_RT_RN_projections(SV2_rel_pos, SV3_rel_pos, 'Non-linear EOMs relative orbits', 'figures/PS3/RTN_projections_numerical.png');

    %plot_RT_RN_projections_both(SV2_HCW_pos,SV2_rel_pos,SV3_HCW_pos,SV3_rel_pos, 'Comparing HCW with non-linear propagation', 'figures/PS3/RTN_projections_HCW_comparison.png');
    %plot_3D_rel_orbit(SV2_HCW_pos,SV3_HCW_pos, 'figures/PS3/3D_HCW_orbit_SV2.png', 'figures/PS3/3D_HCW_orbit_SV3.png');

    % %%%%% PLOT OF RELATIVE ORBITS USING YA %%%%%
    % fprintf('Average rho-to-r_o ratio for SV2: %.3e\n',rho_pos_ratio_SV2_avg);
    % fprintf('Average rho-to-r_o ratio for SV3: %.3e\n',rho_pos_ratio_SV3_avg);
    % fprintf('K_YA for SV2: \n');
    % disp(K_YA_SV2);
    % fprintf('K_YA for SV3: \n');
    % disp(K_YA_SV3);
    % fprintf('Unscaled initial ROE for SV2: \n');
    % disp(ROE_SV2_unscaled);
    % fprintf('Unscaled initial ROE for SV3: \n');
    % disp(ROE_SV3_unscaled); 
    % 
    % % Non-linear (ground truth)
    % plot_rel_pos_vel_single(t_2, t_orbit, SV2_rel_pos_ecc, SV2_rel_vel_ecc, 'Non-linear method: SV2 relative position + velocity with eccentric chief', 'figures/PS3/ecc_nonlinear_pos_vel_SV2.png');
    % plot_rel_pos_vel_single(t_2, t_orbit, SV3_rel_pos_ecc, SV3_rel_vel_ecc, 'Non-linear method: SV3 relative position + velocity with eccentric chief', 'figures/PS3/ecc_nonlinear_pos_vel_SV3.png');
    % plot_RT_RN_projections(SV2_rel_pos_ecc, SV3_rel_pos_ecc, 'Non-linear method: relative orbits with eccentric chief', 'figures/PS3/RTN_projections_numerical_eccentric.png');
    % plot_3D_rel_orbit(SV2_rel_vel_ecc,SV2_rel_vel_ecc,'Non-linear method: 3D orbits','figures/PS3/ecc_nonlinear_3D_orbits.png')
    % 
    % % YA Solution
    % plot_rel_pos_vel_single(t_2, t_orbit, SV2_YA_pos, SV2_YA_vel, 'YA solution: SV2 relative position + velocity', 'figures/PS3/YA_pos_vel_SV2.png');
    % plot_rel_pos_vel_single(t_2, t_orbit, SV3_YA_pos, SV3_YA_vel, 'YA solution: SV3 relative position + velocity', 'figures/PS3/YA_pos_vel_SV3.png');
    % plot_RT_RN_projections(SV2_YA_pos, SV3_YA_pos, 'YA solution: relative orbits', 'figures/PS3/RTN_projections_YA.png');
    % plot_3D_rel_orbit(SV2_YA_pos,SV2_YA_vel,'YA solution: 3D orbits','figures/PS3/YA_sol_3D_orbits.png')
    % 
    % % Geometric mapping
    % plot_rel_pos_vel_single(t_2, t_orbit, SV2_YA_mapping_pos, SV2_YA_mapping_vel, 'Geometric mapping: SV2 relative position + velocity', 'figures/PS3/YA_mapping_pos_vel_SV2.png');
    % plot_rel_pos_vel_single(t_2, t_orbit, SV3_YA_mapping_pos, SV3_YA_mapping_vel, 'Geometric mapping: SV3 relative position + velocity', 'figures/PS3/YA_mapping_pos_vel_SV3.png');
    % plot_RT_RN_projections(SV2_YA_mapping_pos, SV3_YA_mapping_pos, 'Geometric mapping: relative orbits', 'figures/PS3/RTN_projections_YA_mapping.png');
    % plot_3D_rel_orbit(SV2_YA_mapping_pos,SV2_YA_mapping_vel,'Geometric mapping: 3D orbits','figures/PS3/YA_mapping_3D_orbits.png')
    % 
    % % Combined plots
    plot_RT_RN_projections_triple(SV2_YA_pos,SV2_YA_mapping_pos,SV2_rel_pos_ecc,SV3_YA_pos,SV3_YA_mapping_pos,SV3_rel_pos_ecc, 'Comparison of YA methods', 'figures/PS3/RTN_projections_YA_comparison.png');
    plot_RT_RN_error_projections(SV2_YA_pos,SV2_YA_mapping_pos,SV2_rel_pos_ecc,SV3_YA_pos,SV3_YA_mapping_pos,SV3_rel_pos_ecc,'Comparison of errors in YA methods','figures/PS3/RTN_error_projections_YA_comparison.png')
    plot_3D_rel_orbits_triple(SV2_YA_pos,SV2_YA_mapping_pos,SV2_rel_pos_ecc,SV3_YA_pos,SV3_YA_mapping_pos,SV3_rel_pos_ecc,'Comparison of 3D orbits','figures/PS3/3D_YA_comparison.png');

    % % Combined plots with non-zero semi-major axis separation
    % plot_RT_RN_projections_triple(SV2_YA_pos_2,SV2_YA_mapping_pos_2,SV2_rel_pos_ecc_2,SV3_YA_pos_2,SV3_YA_mapping_pos_2,SV3_rel_pos_ecc_2, 'Comparison of YA methods with semi-major axis separation', 'figures/PS3/RTN_projections_YA_comparison_2.png');
    % plot_RT_RN_error_projections(SV2_YA_pos_2,SV2_YA_mapping_pos_2,SV2_rel_pos_ecc_2,SV3_YA_pos_2,SV3_YA_mapping_pos_2,SV3_rel_pos_ecc_2,'Comparison of errors in YA methods with semi-major axis separation','figures/PS3/RTN_error_projections_YA_comparison_2.png')
    % plot_3D_rel_orbits_triple(SV2_YA_pos_2,SV2_YA_mapping_pos_2,SV2_rel_pos_ecc_2,SV3_YA_pos_2,SV3_YA_mapping_pos_2,SV3_rel_pos_ecc_2,'Comparison of 3D orbits with semi-major axis separation','figures/PS3/3D_YA_comparison_2.png');
    % 
    % % Combined plots with highly-eccentric chief orbit
    % plot_RT_RN_projections_triple(SV2_YA_pos_3,SV2_YA_mapping_pos_3,SV2_rel_pos_ecc_3,SV3_YA_pos_3,SV3_YA_mapping_pos_3,SV3_rel_pos_ecc_3, 'Comparison of YA methods with highly-eccentric chief orbit', 'figures/PS3/RTN_projections_YA_comparison_3.png');
    % plot_RT_RN_error_projections(SV2_YA_pos_3,SV2_YA_mapping_pos_3,SV2_rel_pos_ecc_3,SV3_YA_pos_3,SV3_YA_mapping_pos_3,SV3_rel_pos_ecc_3,'Comparison of errors in YA methods with highly-eccentric chief orbit','figures/PS3/RTN_error_projections_YA_comparison_3.png')
    % plot_3D_rel_orbits_triple(SV2_YA_pos_3,SV2_YA_mapping_pos_3,SV2_rel_pos_ecc_3,SV3_YA_pos_3,SV3_YA_mapping_pos_3,SV3_rel_pos_ecc_3,'Comparison of 3D orbits with highly-eccentric chief orbit','figures/PS3/3D_YA_comparison_3.png');

    %plot_rel_pos_vel_single(t_2, t_orbit, SV2_YA_mapping_pos, SV2_YA_mapping_vel, 'SV2 YA gemetric mapping relative position + velocity', 'figures/PS3/YA_mapping_pos_vel_SV2.png');
    %plot_rel_pos_vel_single(t_2, t_orbit, SV3_YA_mapping_pos, SV3_YA_mapping_vel, 'SV3 YA gemetric mapping relative position + velocity', 'figures/PS3/YA_mapping_pos_vel_SV3.png');
    %plot_RT_RN_projections(SV2_YA_mapping_pos, SV3_YA_mapping_pos, 'YA geometric mapping relative orbits', 'figures/PS3/RTN_projections_YA_mapping.png');
   
    %plot_3D_rel_orbits_both(SV2_YA_pos,SV2_YA_mapping_pos,SV3_YA_pos,SV3_YA_mapping_pos,'figures/PS3/3D_YA_comparison.png');

end