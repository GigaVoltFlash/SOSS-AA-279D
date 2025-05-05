%% AA279D Simulations for PSET 5
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

%%%%% INITIAL RELATIVE ORBITAL ELEMENTS %%%%%%% 
% From PS3
d_a_SV2_init_1 = 0; % m
d_lambda_SV2_init_1 = 0; % m % -124
d_e_x_SV2_init_1 = 0; % m
d_e_y_SV2_init_1 = 100; % m
d_i_x_SV2_init_1 = 0; % m % 79
d_i_y_SV2_init_1 = 1000; % m % 1005

d_a_SV3_init_1 = 0; % m
d_lambda_SV3_init_1 = 0; % m % -793
d_e_x_SV3_init_1 = 0; % m
d_e_y_SV3_init_1 = 200; % m
d_i_x_SV3_init_1 = 0; % m
d_i_y_SV3_init_1 = 800; % m % 827

plot_rmin_contour(a_SV1_init*1e3, [10, 250, 1000]);

% SIM IN RTN?
% Need a high fidelity sim with J2 perturbations that we can input
% different control modes into.
% Features of this simulation
%       * Take in functions for the different modes that define control.
%       * Apply specific delta v maneuvers at times specified
%       * Output RTN space throughout simulation but also in ROE space (or
%       maybe do conversion post)
%       * Be able to simulate 1 to 4 modes (as in turn, off particular
%       modes)
%       * Intelligently activate the different modes based on a state
%       machine (i.e. or just hard-code when we want specific maneuvers).
