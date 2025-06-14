%% AA279D Simulations for PSET 7
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

% Common initialization file for full sim with all modes
initialize_roe_and_modes;

test_EKF_continuous(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
    SV1_OE_init, state_abs_SV2_init, state_abs_SV3_init , ...
    state_init, t_orbit, t_series);