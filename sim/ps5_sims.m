%% AA279D Simulations for PSET 5
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

%%%%% DEFINE ROE FOR EACH MODE %%%%%%% 
%plot_rmin_contour(a_SV1_init*1e3, [10, 250, 1000]);



% d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y
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


for mode = 1:4
    SV2_ROE = SV2_modes(mode, :);
    SV3_ROE = SV3_modes(mode, :);

    % Initial absolute ECI states of deputies
    [r_SV2_init, v_SV2_init]  = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
        SV2_ROE(1), SV2_ROE(2), SV2_ROE(3), SV2_ROE(4), SV2_ROE(5), SV2_ROE(6));
    [r_SV3_init, v_SV3_init] = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
        SV3_ROE(1), SV3_ROE(2), SV3_ROE(3), SV3_ROE(4), SV3_ROE(5), SV3_ROE(6));

    state_abs_SV2_init = [r_SV2_init; v_SV2_init];
    state_abs_SV3_init = [r_SV3_init; v_SV3_init];

    fig_path = sprintf('figures/PS5/mode_%d_RTN.png', mode);
    title_str = sprintf('ROE Mode %d: SV2 & SV3 Relative Motion', mode);

    simulate_and_plot_relative_motion(t_series, SV2_ROE, state_abs_SV2_init, ...
                                      SV3_ROE, state_abs_SV3_init, SV1_OE_init, ...
                                      state_init, fig_path, title_str);
end