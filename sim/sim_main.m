%% AA279D Master Simulation
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 
close all;
addpath("plotting_methods\");
addpath("mean_osc\");

% Load all constants from the constants file
constants()
global mu_earth

a_init = 6944;
ex_init = -4e-5;
ey_init = 1.6e-3;
i_init = 99.4;
RAAN_init = -151.1;
u_init = -47.9;
[a_init, e_init, i_init, RAAN_init, w_init, nu_init, M_init] = quasi_nonsing2OE(a_init, ex_init, ey_init, i_init, RAAN_init, u_init);

[r_init, v_init] = OE2ECI(a_init, e_init, i_init, RAAN_init, w_init, nu_init);
state_init = [r_init; v_init];
tstart = 0; % Start time (s)

t_orbit = 2*pi*sqrt(a_init^3/mu_earth);
tend = 25*t_orbit; % End time (s)
tint = t_orbit/50.0; % Time step (s)

%%%%%%%%%%% RUN SIM WITH AND WITHOUT J2 %%%%%%%%%%%
[t_1, state1] = rk4_eom_ECI(tstart:tint:tend, state_init, false);
[t_2, state2] = rk4_eom_ECI(tstart:tint:tend, state_init, true);
r_ECI_with_j2 = state2(:,1:3);
v_ECI_with_j2 = state2(:,4:6);
r_ECI_no_j2 = state1(:,1:3);
v_ECI_no_j2 = state1(:,4:6);

%%%%%% RUN SIM OF RELATIVE MOTION %%%%%%%%
[t_3, state3] = rk4_eom_rel_RTN(tstart:tint:tend, state_rel_SV2_init);
[t_4, state4] = rk4_eom_rel_RTN(tstart:tint:tend, state_rel_SV3_init);
r_ECI_no_j2_2 = state3(:,1:3); % This is to check that the propagation gives the same results as the previous propagation
v_ECI_no_j2_2 = state3(:,4:6);
SV2_rel_pos = [state3(:, 7), state3(:, 8), state3(:, 9)];
SV3_rel_pos = [state4(:, 7), state4(:, 8), state4(:, 9)];

%%%%%%% KEPLERIAN PROPAGATION %%%%%%% 
% Propagates over the same times t_1 as the sim without J2 
[r_ECI_keplerian, v_ECI_keplerian] = keplerian_propagator(a_init, e_init, i_init, RAAN_init, w_init, t_1, M_init);

%%%%%% CONVERT TO RTN %%%%%%%%
[r_RTN_no_j2, v_RTN_no_j2] = ECI2RTN(r_ECI_no_j2,v_ECI_no_j2);
[r_RTN_keplerian, v_RTN_keplerian] = ECI2RTN(r_ECI_keplerian,v_ECI_keplerian);

%%%%% INITIAL RELATIVE ORBITAL ELEMENTS %%%%%%%
d_a_2_init = 21; % m
d_lambda_2_init = -124350; % m
d_e_x_2_init = 110; % m
d_e_y_2_init = 202; % m
d_i_x_2_init = 79; % m
d_i_y_2_init = 1005; % m

[a_2_init,e_x_2_init,e_y_2_init,i_2_init,RAAN_2_init,u_2_init] = ...
ROE2quasi_nonsing(a_init,ex_init,ey_init,i_init,RAAN_init,u_init, ...
d_a_2_init,d_lambda_2_init,d_e_x_2_init,d_e_y_2_init,d_i_x_2_init,d_i_y_2_init);

[a_2_init, e_2_init, i_2_init, RAAN_2_init, w_2_init, nu_2_init, M_2_init] = ...
quasi_nonsing2OE(a_2_init, e_x_2_init, e_y_2_init, i_2_init, RAAN_2_init, u_2_init);

[r_2_init, v_2_init] = OE2ECI(a_2_init, e_2_init, i_2_init, RAAN_2_init, w_2_init, nu_2_init);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
%%%%%% PLOT 3D COMPARISON PLOTS %%%%%% 
plot_3D_orbit_aoe_compare(r_ECI_no_j2, r_ECI_with_j2, 'Orbit without J2', 'Orbit with J2')
plot_3D_orbit_aoe_compare(r_ECI_no_j2, r_ECI_keplerian, 'Orbit without J2', 'Keplerian Propagator')

%%%%%% COMPARE RTN VALUES %%%%%% 
plot_rtn_compare(t_1, r_RTN_no_j2, r_RTN_keplerian, v_RTN_no_j2, v_RTN_keplerian, 'Numerical with no J2', 'Keplerian', 'Error between No-J2 Numerical and Keplerian, large time-step');

%%%%%%% COMPUTE AND PLOT OE, ECCENTRICITY, ANGULAR MOMENTUM, SPECIFIC MECHANICAL ENERGY %%%%%%% 
[a_2,e_2,i_2,RAAN_2,omega_2,nu_2] = compute_and_plot_orbit_params(r_ECI_no_j2,v_ECI_no_j2,r_ECI_with_j2,v_ECI_with_j2,t_1,t_2);

%%%%%%%% COMPUTE AND PLOT OSCULATING VERSUS MEAN %%%%%%%%%%%%
[a_3,e_3,i_3,RAAN_3,omega_3,nu_3,t_3] = compute_mean_oe(a_init,ex_init,ey_init,i_init,RAAN_init,u_init,tstart,tint,tend);
plot_osc_mean_oe(a_2,e_2,i_2,RAAN_2,omega_2,nu_2,t_2,a_3,e_3,i_3,RAAN_3,omega_3,nu_3,t_3);

