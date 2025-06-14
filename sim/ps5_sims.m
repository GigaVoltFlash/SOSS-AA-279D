%% AA279D Simulations for PSET 5
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

%%%%% DEFINE ROE FOR EACH MODE %%%%%%% 
plot_rmin_contour(a_SV1_init*1e3, [10, 250, 1000]);

% Common initialization file for full sim with all modes
initialize_roe_and_modes;

% for mode = 1:3
%     SV2_ROE = SV2_modes(mode, :);
%     SV3_ROE = SV3_modes(mode, :);
% 
%     % Initial absolute ECI states of deputies
%     [r_SV2_init, v_SV2_init]  = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
%         SV2_ROE(1), SV2_ROE(2), SV2_ROE(3), SV2_ROE(4), SV2_ROE(5), SV2_ROE(6));
%     [r_SV3_init, v_SV3_init] = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init,...
%         SV3_ROE(1), SV3_ROE(2), SV3_ROE(3), SV3_ROE(4), SV3_ROE(5), SV3_ROE(6));
% 
%     state_abs_SV2_init = [r_SV2_init; v_SV2_init];
%     state_abs_SV3_init = [r_SV3_init; v_SV3_init];
% 
%     fig_path = sprintf('figures/PS5/mode_%d_RTN.png', mode);
%     title_str = sprintf('Mode %d: SV2 & SV3 Relative Motion in RTN', mode);
%     fig_path2 = sprintf('figures/PS5/mode_%d_ROE_Planes.png', mode);
%     title_str2 = sprintf('Mode %d: SV2 & SV3 Relative Motion in ROE Planes', mode);
%     fig_path3 = sprintf('figures/PS5/mode_%d_ROE_Time.png', mode);
%     title_str3 = sprintf('Mode %d: SV2 & SV3 Relative Motion in ROE Time History', mode);
% 
%     simulate_and_plot_relative_motion(t_series, t_orbit, SV2_ROE, state_abs_SV2_init, ...
%                                       SV3_ROE, state_abs_SV3_init, SV1_OE_init, ...
%                                       state_init, fig_path, title_str, fig_path2, title_str2, fig_path3, title_str3);
% 
%     roe_initial = SV3_modes(mode,:)/a_SV1_init; % must be unscaled!
%     roe_final = SV3_modes(mode+1,:)/a_SV1_init;
%     %roe_initial = SV3_modes(2,:)/a_SV1_init; % must be unscaled!
%     %roe_final = SV3_modes(3,:)/a_SV1_init;
%     init_time = tstart;
%     final_time = 2*t_orbit;
% 
%     [delta_v_vals, delta_v_times] = mode2_control(roe_initial, roe_final, init_time, final_time, SV1_OE_init, u_SV1_init);
% 
%     fig_path = sprintf('figures/PS5/mode_%d_RTN_maneuvers.png', mode);
%     title_str = sprintf('ROE Mode %d: SV2 & SV3 Relative Motion with manuever', mode);
% end

%%%%%%%%%%%% INTEGRATED SIM %%%%%%%%%%%%%%%
sim_all_maneuvers_sk_impulsive(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
                                     SV1_OE_init, state_abs_SV2_init, state_abs_SV3_init , ...
                                     state_init, t_orbit, t_series,'figures/PS5/RTN_3d_projections_all_maneuvers.png', '');
