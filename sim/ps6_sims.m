%% AA279D Simulations for PSET 6
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

%%%%% DEFINE ROE FOR EACH MODE %%%%%%% 
% (d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y)*a_c in meters
SV2_modes = [
    0, 0, 0,   300, 0, 300;
    0, 0, 0,   300, 0, 300;
    0, 0, 0,   300, 0, 300;
    0, 0, 0,   300, 0, 300
];

SV3_modes = [
    0, 0, 0,   250, 0, -250;
    0, 0, 0,   100, 0, -100;
    0, 0, 0,   10, 0, -10;
    0, 0, 0,    1, 0, -1
];

t_series = tstart:tint:tend;
SV1_OE_init = [a_SV1_init, e_SV1_init, i_SV1_init, RAAN_SV1_init, w_SV1_init, M_SV1_init];
a_chief_meters = a_SV1_init*1e3;

