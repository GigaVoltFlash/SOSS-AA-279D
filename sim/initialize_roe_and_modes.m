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

% Add up to 25 orbits
num_orbits_modes = [0,0,0,0];
num_orbits_station_keep = [5,0,0,0];

SV2_ROE = SV2_modes(1, :);
SV3_ROE = SV3_modes(1, :);

% Initial absolute ECI states of deputies
[r_SV2_init, v_SV2_init]  = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
SV2_ROE(1), SV2_ROE(2), SV2_ROE(3), SV2_ROE(4), SV2_ROE(5), SV2_ROE(6));
[r_SV3_init, v_SV3_init] = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
SV3_ROE(1), SV3_ROE(2), SV3_ROE(3), SV3_ROE(4), SV3_ROE(5), SV3_ROE(6));

state_abs_SV2_init = [r_SV2_init; v_SV2_init];
state_abs_SV3_init = [r_SV3_init; v_SV3_init];

J2_flag = true;

[d_a_SV3_osc, d_lambda_SV3_osc, d_e_x_SV3_osc, d_e_y_SV3_osc, d_i_x_SV3_osc, d_i_y_SV3_osc, ...
         d_a_SV3_mean, d_lambda_SV3_mean, d_e_x_SV3_mean, d_e_y_SV3_mean, d_i_x_SV3_mean, d_i_y_SV3_mean, ...
         a_osc, e_x_osc, e_y_osc, i_osc, RAAN_osc, u_osc, ...
         a_mean, e_x_mean, e_y_mean, i_mean, RAAN_mean, u_mean] = ...
            compute_OE_ROE_mean_osc(r_SV1_ECI_init', v_SV1_ECI_init', r_SV3_init', v_SV3_init', J2_flag);

[r_SV3_mean_init, v_SV3_mean_init] = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
d_a_SV3_osc, d_lambda_SV3_osc, d_e_x_SV3_osc, d_e_y_SV3_osc, d_i_x_SV3_osc, d_i_y_SV3_osc);

state_abs_SV3_mean_init = [r_SV3_mean_init; v_SV3_mean_init];