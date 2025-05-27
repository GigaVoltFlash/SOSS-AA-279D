%% AA279D Simulations for PSET 8
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

% Common initialization file for full sim with all modes
initialize_roe_and_modes;

% using SV3 mean abs state init
EKF_continuous_no_control(SV1_OE_init, state_abs_SV2_init, state_abs_SV3_init , state_init, t_orbit, t_series);