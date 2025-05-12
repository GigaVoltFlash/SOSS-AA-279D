%% AA279D Simulations for PSET 5
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

%%%%% DEFINE ROE FOR EACH MODE %%%%%%% 
plot_rmin_contour(a_SV1_init*1e3, [10, 250, 1000]);



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

for mode = 1:2
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
    title_str = sprintf('Mode %d: SV2 & SV3 Relative Motion in RTN', mode);
    fig_path2 = sprintf('figures/PS5/mode_%d_ROE_Planes.png', mode);
    title_str2 = sprintf('Mode %d: SV2 & SV3 Relative Motion in ROE Planes', mode);
    fig_path3 = sprintf('figures/PS5/mode_%d_ROE_Time.png', mode);
    title_str3 = sprintf('Mode %d: SV2 & SV3 Relative Motion in ROE Time History', mode);

    simulate_and_plot_relative_motion(t_series, t_orbit, SV2_ROE, state_abs_SV2_init, ...
                                      SV3_ROE, state_abs_SV3_init, SV1_OE_init, ...
                                      state_init, fig_path, title_str, fig_path2, title_str2, fig_path3, title_str3);

    roe_initial = SV3_modes(mode,:)/a_SV1_init; % must be unscaled!
    roe_final = SV3_modes(mode+1,:)/a_SV1_init;
    %roe_initial = SV3_modes(2,:)/a_SV1_init; % must be unscaled!
    %roe_final = SV3_modes(3,:)/a_SV1_init;
    init_time = tstart;
    final_time = 2*t_orbit;

    [delta_v_vals, delta_v_times] = mode2_control(roe_initial, roe_final, init_time, final_time, SV1_OE_init, u_SV1_init);

    fig_path = sprintf('figures/PS5/mode_%d_RTN_maneuvers.png', mode);
    title_str = sprintf('ROE Mode %d: SV2 & SV3 Relative Motion with manuever', mode);

    %simulate_and_plot_relative_motion_with_maneuvers(t_series, t_orbit, SV2_ROE, state_abs_SV2_init, ...
    %    SV3_ROE, state_abs_SV3_init, SV1_OE_init, state_init, ...
    %    fig_path, title_str, delta_v_times, delta_v_vals)
end

num_orbits_modes = [0,3,3,3];
num_orbits_station_keep = [5,4,4,3];


SV2_ROE = SV2_modes(1, :);
SV3_ROE = SV3_modes(1, :);

% Initial absolute ECI states of deputies
[r_SV2_init, v_SV2_init]  = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
SV2_ROE(1), SV2_ROE(2), SV2_ROE(3), SV2_ROE(4), SV2_ROE(5), SV2_ROE(6));
[r_SV3_init, v_SV3_init] = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
SV3_ROE(1), SV3_ROE(2), SV3_ROE(3), SV3_ROE(4), SV3_ROE(5), SV3_ROE(6));

state_abs_SV2_init = [r_SV2_init; v_SV2_init];
state_abs_SV3_init = [r_SV3_init; v_SV3_init];

%%%%%%%%% DEPRECATED CODE, NEW INTEGRATED SIM BELOW %%%%%%%%%%
% sim_and_plot_control_modes_old(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
%                                      SV1_OE_init, u_SV1_init, state_abs_SV2_init, state_abs_SV3_init , ...
%                                      state_init, t_orbit, t_series, tint);
% sim_and_plot_control_modes(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
%                                      SV1_OE_init, u_SV1_init, state_abs_SV2_init, state_abs_SV3_init , ...
%                                      state_init, t_orbit, t_series, tint, 'figures/PS5/control_mode_test.png', 'Here we go');

% % Starting at mode 2 for the station keeping system
% SV2_ROE = SV2_modes(2, :);
% SV3_ROE = SV3_modes(2, :);
% 
% % Initial absolute ECI states of deputies
% [r_SV2_init, v_SV2_init]  = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
% SV2_ROE(1), SV2_ROE(2), SV2_ROE(3), SV2_ROE(4), SV2_ROE(5), SV2_ROE(6));
% [r_SV3_init, v_SV3_init] = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
% SV3_ROE(1), SV3_ROE(2), SV3_ROE(3), SV3_ROE(4), SV3_ROE(5), SV3_ROE(6));
% 
% state_abs_SV2_init = [r_SV2_init; v_SV2_init];
% state_abs_SV3_init = [r_SV3_init; v_SV3_init];
% 
% sim_and_plot_stationkeeping(SV2_modes, SV3_modes, SV1_OE_init, state_abs_SV2_init, state_abs_SV3_init, ...
%                                      state_init, t_orbit, t_series, 'figures/PS5/station_keeping.png', 'Testing station keeping');

sim_all_maneuvers_station_keeping(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
                                     SV1_OE_init, state_abs_SV2_init, state_abs_SV3_init , ...
                                     state_init, t_orbit, t_series,'figures/PS5/RTN_3d_projections_all_maneuvers.png', '');
