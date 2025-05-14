function sim_all_maneuvers_sk_continuous(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, SV1_OE_init, SV2_state_init, SV3_state_init, ...
                                     SV1_state_init, t_orbit, t_series, fig_path, title_str)
    full_times = t_series(:);
    a_chief = SV1_OE_init(1);

    % Factors for Lyaponuv feedback gain matrix  
    N = 14;
    k = 1000;

    switch_times = t_orbit*cumsum(reshape([num_orbits_modes; num_orbits_station_keep], 8, 1));

    % Initialize with same size as full time vector
    state_SV1_all = zeros(length(t_series), 6);
    state_SV2_all = zeros(length(t_series), 6);
    state_SV3_all = zeros(length(t_series), 6);

    % Nominal ROE for particular modes
    SV2_roe_nom = SV2_modes(2, :)/a_chief;
    SV3_roe_nom_mode1 = SV3_modes(1, :)/a_chief;
    SV3_roe_nom_mode2 = SV3_modes(2, :)/a_chief;
    SV3_roe_nom_mode3 = SV3_modes(3, :)/a_chief;
    SV3_roe_nom_mode4 = SV3_modes(4, :)/a_chief;
    SV3_roe_nom = SV3_roe_nom_mode1;

    % Maybe pick these differently later
    SV2_delta_di_max = 0.5/a_chief;
    SV3_delta_di_max = 0.5/a_chief;
    SV2_delta_de_max = 0.5/a_chief;
    SV3_delta_de_max = 0.5/a_chief;

    % Setting initial states
    SV1_state = SV1_state_init;
    SV2_state = SV2_state_init;
    SV3_state = SV3_state_init;

    state_SV3_all(1, :) = SV3_state';
    state_SV2_all(1, :) = SV2_state';
    state_SV1_all(1, :) = SV1_state';

    dt = full_times(2) - full_times(1);
    n_steps = length(full_times);

    % Initialize delta-v and acceleration storage
    SV2_dv_vals = zeros(length(t_series), 3);
    SV3_dv_vals = zeros(length(t_series), 3);
    SV2_a_vals = zeros(length(t_series), 3);
    SV3_a_vals = zeros(length(t_series), 3);

    withJ2 = true;

    for i=2:n_steps
        t = full_times(i-1);

        % PRIMARY PROPAGATION (RK4 STEPS)
        k1_SV1 = eom_ECI(t, SV1_state, withJ2);
        k2_SV1 = eom_ECI(t + dt/2, SV1_state + dt/2 * k1_SV1, withJ2);
        k3_SV1 = eom_ECI(t + dt/2, SV1_state + dt/2 * k2_SV1, withJ2);
        k4_SV1 = eom_ECI(t + dt,   SV1_state + dt   * k3_SV1, withJ2);
        SV1_state = SV1_state + dt/6 * (k1_SV1 + 2*k2_SV1 + 2*k3_SV1 + k4_SV1);
        
        k1_SV2 = eom_ECI(t, SV2_state, withJ2);
        k2_SV2 = eom_ECI(t + dt/2, SV2_state + dt/2 * k1_SV2, withJ2);
        k3_SV2 = eom_ECI(t + dt/2, SV2_state + dt/2 * k2_SV2, withJ2);
        k4_SV2 = eom_ECI(t + dt,   SV2_state + dt   * k3_SV2, withJ2);
        SV2_state = SV2_state + dt/6 * (k1_SV2 + 2*k2_SV2 + 2*k3_SV2 + k4_SV2);
        
        k1_SV3 = eom_ECI(t, SV3_state, withJ2);
        k2_SV3 = eom_ECI(t + dt/2, SV3_state + dt/2 * k1_SV3, withJ2);
        k3_SV3 = eom_ECI(t + dt/2, SV3_state + dt/2 * k2_SV3, withJ2);
        k4_SV3 = eom_ECI(t + dt,   SV3_state + dt   * k3_SV3, withJ2);
        SV3_state = SV3_state + dt/6 * (k1_SV3 + 2*k2_SV3 + 2*k3_SV3 + k4_SV3);

        % Save states to buffer
        state_SV3_all(i, :) = SV3_state';
        state_SV2_all(i, :) = SV2_state';
        state_SV1_all(i, :) = SV1_state';

        % Convert current states to ROE
        [d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2] = ECI2ROE_array_mean(SV1_state(1:3)', SV1_state(4:6)', SV2_state(1:3)', SV2_state(4:6)', true);
        [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3] = ECI2ROE_array_mean(SV1_state(1:3)', SV1_state(4:6)', SV3_state(1:3)', SV3_state(4:6)', true);

        [a,e,inc,RAAN,omega,nu,M] = ECI2OE(SV1_state(1:3), SV1_state(4:6));
        SV1_oe = [a,e,inc,RAAN,omega,nu,M];
        SV2_roe = [d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2]/a_chief;
        SV3_roe = [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3]/a_chief;

        % SV2 control inputs are always determined the same way
        SV2_a_vals(i, :) = station_keeping_continuous(SV2_roe, SV2_roe_nom, SV2_delta_de_max, SV2_delta_di_max, SV1_oe, N, k, 1/t_orbit*switch_times(end));

        % This is primarily for SV3, SV2 just continuously station keeps     
        if t < switch_times(1)
            % do nothing, already at initial situation
        elseif t < switch_times(3) && t > switch_times(2)
            SV3_roe_nom = SV3_roe_nom_mode2; % Switches it for maneuver and station keeping
            SV3_a_vals(i, :) = Lyapunov_feedback_control(SV3_roe, SV3_roe_nom, SV1_oe, N, k, num_orbits_modes(2));
        elseif t < switch_times(5) && t > switch_times(4)
            SV3_roe_nom = SV3_roe_nom_mode3; % Switches it for maneuver and station keeping
            SV3_a_vals(i, :) = Lyapunov_feedback_control(SV3_roe, SV3_roe_nom, SV1_oe, N, k, num_orbits_modes(2));
        elseif t < switch_times(7) && t > switch_times(6)
            SV3_roe_nom = SV3_roe_nom_mode4; % Switches it for maneuver and station keeping
            SV3_a_vals(i, :) = Lyapunov_feedback_control(SV3_roe, SV3_roe_nom, SV1_oe, N, k, num_orbits_modes(2));
        else
            % Run station keeping 
            % SV3_a_vals(i, :) = station_keeping_continuous(SV3_roe, SV3_roe_nom, SV3_delta_de_max, SV3_delta_di_max, SV1_oe, N, k);
            SV3_a_vals(i, :) = Lyapunov_feedback_control(SV3_roe, SV3_roe_nom, SV1_oe, N, k, num_orbits_station_keep(2)); % TODO: Fix this hardcoding
            % SV3_a_vals(i, :) = [0, 0, 0];
        end

        %%% APPLY ACCELERATION AS DELTA V 
        % Convert a_RTN to dv_RTN by multiplying by dt
        SV2_dv_vals(i,:) = dv_RTN2ECI(SV1_state(1:3), SV1_state(4:6), dt*SV2_a_vals(i, :)'/1e3); % m/s^2 --> m/s --> km/s
        SV3_dv_vals(i,:) = dv_RTN2ECI(SV1_state(1:3), SV1_state(4:6), dt*SV3_a_vals(i, :)'/1e3); % m/s^2 --> m/s --> km/s

        SV2_state(4:6) = SV2_state(4:6) + SV2_dv_vals(i,:)';
        SV3_state(4:6) = SV3_state(4:6) + SV3_dv_vals(i,:)';

    end

    % Final SV1 propagated states
    r_SV1 = state_SV1_all(:,1:3);
    v_SV1 = state_SV1_all(:,4:6);
    % Final SV2 propagated states
    r_SV2 = state_SV2_all(:,1:3);
    v_SV2 = state_SV2_all(:,4:6);
    % Final SV3 propagated states
    r_SV3 = state_SV3_all(:,1:3);
    v_SV3 = state_SV3_all(:,4:6);

    % --- Relative RTN position computation ---
    SV2_rel_pos = zeros(length(t_series), 3);
    SV3_rel_pos = zeros(length(t_series), 3);

    for i = 1:length(t_series)
        [rho2, ~] = ECI2RTN_rel(r_SV1(i,:), v_SV1(i,:), r_SV2(i,:), v_SV2(i,:));
        [rho3, ~] = ECI2RTN_rel(r_SV1(i,:), v_SV1(i,:), r_SV3(i,:), v_SV3(i,:));
        SV2_rel_pos(i,:) = rho2;
        SV3_rel_pos(i,:) = rho3;
    end

    % --- Plot ---
    plot_3D_rel_orbit(SV2_rel_pos, SV3_rel_pos, true, 'All maneuvers 3D path', 'figures/PS6/ROE_3d_all_maneuvers.png');
    % create_roe_3d_animation(SV2_rel_pos, SV3_rel_pos, 'SV2 and SV3 3D Relative Orbits (Animated)');
    
    RT = true; RN = true; NT = true;
    plot_RT_RN_projections_separate(SV2_rel_pos, SV3_rel_pos, RT, RN, NT, title_str, fig_path);
    [d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2] = ECI2ROE_array_mean(r_SV1, v_SV1, r_SV2, v_SV2, true);
    [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3] = ECI2ROE_array_mean(r_SV1, v_SV1, r_SV3, v_SV3, true);
    %plot_ROE_planes(full_times, t_orbit, d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3, 'figures/PS6/ROE_planes_SV3.png', 'figures/PS6/ROE_over_time_SV3.png');
    %plot_ROE_planes(full_times, t_orbit, d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2, 'figures/PS6/ROE_planes_SV2.png', 'figures/PS6/ROE_over_time_SV2.png');

    plot_ROE_planes_with_modes(full_times, t_orbit, d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2,SV2_modes,...
       num_orbits_modes, num_orbits_station_keep,'figures/PS6/ROE_planes_modes_SV2.png', 'figures/PS6/ROE_over_time_modes_SV2.png', 'figures/PS6/ROE_error_over_time_modes_SV2.png');
    plot_ROE_planes_with_modes(full_times, t_orbit, d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3,SV3_modes,...
        num_orbits_modes, num_orbits_station_keep,'figures/PS6/ROE_planes_modes_SV3.png', 'figures/PS6/ROE_over_time_modes_SV3.png', 'figures/PS6/ROE_error_over_time_modes_SV3.png');
    plot_delta_v_timeline(full_times, SV2_a_vals, t_orbit, SV2_modes, num_orbits_modes, num_orbits_station_keep, ...
       true, 'figures/PS6/delta_v_timeline_modes_SV2.png');
    plot_delta_v_timeline(full_times, SV3_a_vals, t_orbit, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
       true, 'figures/PS6/delta_v_timeline_modes_SV3.png');
    plot_delta_v_cumulative_timeline(full_times, SV3_a_vals, t_orbit, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
       true, 'figures/PS6/delta_v_cumulative_timeline_modes_SV3.png')
    
end