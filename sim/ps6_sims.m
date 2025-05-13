%% AA279D Simulations for PSET 6
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

% Common initialization file for full sim with all modes
initialize_roe_and_modes;

sim_all_maneuvers_sk_continuous(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
    SV1_OE_init, state_abs_SV2_init, state_abs_SV3_init , ...
    state_init, t_orbit, t_series,'figures/PS5/RTN_3d_projections_all_maneuvers_cont.png', '');


% test 
SV1_oe = [6.943636694889446e+03,0.001636285638892,99.400248688773600,-1.506974728946212e+02,91.286184249414080,2.195575939795433e+02,-1.403228800441063e+02];
a_chief = SV1_oe(1);
SV3_roe_nom_mode1 = SV3_modes(1, :)/a_chief;
SV3_roe_curr = [9.852476243127120e-05,-0.004729063592322,6.819047991859608e-04,0.035895731935810,-9.701346928956628e-06,-0.035987943224596];
SV3_roe_nom_mode2 = SV3_modes(2, :)/a_chief;

control_vec_result = Lyapunov_feedback_control(SV3_roe_curr, SV3_roe_nom_mode2, SV1_oe, 8, 1e5);
disp(control_vec_result)