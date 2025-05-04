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