%% AA279D Simulations for PSET 9
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

% Common initialization file for full sim with all modes
initialize_roe_and_modes;

EKF_continuous_with_control(SV1_OE_init, SV2_ROE, SV3_ROE, SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, num_orbits_total, t_orbit, t_series);
