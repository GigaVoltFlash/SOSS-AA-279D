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

%plot_rmin_contour(a_SV1_init*1e3, [10, 250, 1000]);

%delta_r_min = 250;         % Desired minimum distance = 50 m
%delta_e = 0.000005;         % Example delta e
%delta_i = 0.000036;

%delta_r_min_calc = sqrt(2)*a_SV1_init*1e3*abs(delta_e*delta_i)/...
%    sqrt(delta_e^2 + delta_i^2 + abs(delta_e+delta_i)*abs(delta_e-delta_i));
%disp(delta_r_min_calc)

%delta_i = solve_delta_i(a_SV1_init*1e3, delta_r_min, delta_e);
%fprintf('|delta i| = %.6f\n', delta_i);

SV2_ROE_Mode_1 = [d_a_SV2_init_1, d_lambda_SV2_init_1, d_e_x_SV2_init_1, d_e_y_SV2_init_1, d_i_x_SV2_init_1, d_i_y_SV2_init_1];
SV3_ROE_Mode_1 = [d_a_SV3_init_1, d_lambda_SV3_init_1, d_e_x_SV3_init_1, d_e_y_SV3_init_1, d_i_x_SV3_init_1, d_i_y_SV3_init_1];
SV1_OE_init = [a_SV1_init, e_SV1_init, i_SV1_init, RAAN_SV1_init, w_SV1_init, M_SV1_init];

t_series = tstart:tint:tend;

simulate_and_plot_relative_motion(t_series, SV2_ROE_Mode_1, SV3_ROE_Mode_1, SV1_OE_init, 'figures/PS5/test_ROE.png', 'Is this working?')