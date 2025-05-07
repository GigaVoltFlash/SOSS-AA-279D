function simulate_and_plot_relative_motion_with_maneuvers(t, t_orbit, SV2_roe_init, state_abs_SV2_init, ...
    SV3_roe_init, state_abs_SV3_init, SV1_oe_init, SV1_state_init, ...
    delta_v_times, delta_v_vals, SV3_modes, num_orbits_modes, num_orbits_station_keep)

    % % Propagate chief (with J2)
    % [~, chief_state] = rk4_eom_ECI(t, SV1_state_init, true);
    % r_SV1 = chief_state(:,1:3);
    % v_SV1 = chief_state(:,4:6);

    % % Propagate SV2 (no impulses)
    % [~, state_SV2] = rk4_eom_ECI(t, state_abs_SV2_init, true);
    % r_SV2 = state_SV2(:,1:3);
    % v_SV2 = state_SV2(:,4:6);

    % --- Split SV3 propagation into segments with impulsive delta-v's ---
    full_times = t(:);
    man_times = [0; delta_v_times(:); full_times(end)];
    % Initialize with same size as full time vector
    state_SV1_all = zeros(length(t), 6);
    state_SV2_all = zeros(length(t), 6);
    state_SV3_all = zeros(length(t), 6);

    
    % Index tracker for filling in state_SV3_all
    state_curr_SV3 = state_abs_SV3_init';
    state_curr_SV2 = state_abs_SV2_init';
    state_curr_SV1 = SV1_state_init';
    
    % Loop over maneuver segments
    % Currently just maneuvering SV3
    for i = 1:length(man_times) - 1
        % Find index range for this time segment
        idx_range = find(t >= man_times(i) & t <= man_times(i+1));
    
        if isempty(idx_range)
            continue; % Skip empty segments (shouldn’t happen with valid times)
        end
    
        % Time vector for this segment
        t_segment = t(idx_range);
    
        % Apply delta-v at start of segment (except for first)
        if i > 1
            dv = dv_RTN2ECI(state_curr_SV1(1:3), state_curr_SV1(4:6), delta_v_vals(i-1, :)'/1e3);  % delta-v at maneuver time (m/s --> km/s)
            state_curr_SV3(4:6) = state_curr_SV3(4:6) + dv;
        end
    
        % Propagate with rk4
        [~, state_segment_SV3] = rk4_eom_ECI(t_segment, state_curr_SV3, true);
        [~, state_segment_SV2] = rk4_eom_ECI(t_segment, state_curr_SV2, true);
        [~, state_segment_SV1] = rk4_eom_ECI(t_segment, state_curr_SV1, true);
    
        % Save to full array
        state_SV3_all(idx_range, :) = state_segment_SV3;
        state_SV2_all(idx_range, :) = state_segment_SV2;
        state_SV1_all(idx_range, :) = state_segment_SV1;
    
        % Update for next segment
        state_curr_SV3 = state_segment_SV3(end, :)';
        state_curr_SV2 = state_segment_SV2(end, :)';
        state_curr_SV1 = state_segment_SV1(end, :)';
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
    SV2_rel_pos = zeros(length(t), 3);
    SV3_rel_pos = zeros(length(t), 3);

    for i = 1:length(t)
        [rho2, ~] = ECI2RTN_rel(r_SV1(i,:), v_SV1(i,:), r_SV2(i,:), v_SV2(i,:));
        [rho3, ~] = ECI2RTN_rel(r_SV1(i,:), v_SV1(i,:), r_SV3(i,:), v_SV3(i,:));
        SV2_rel_pos(i,:) = rho2;
        SV3_rel_pos(i,:) = rho3;
    end

    % --- Plot ---
    RT = true; RN = true; NT = true;
    plot_RT_RN_projections_separate(SV2_rel_pos, SV3_rel_pos, RT, RN, NT, 'RTN Control Modes', 'figures/PS5/RTN_control_modes.png');

    [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3] = ECI2ROE_array_mean(r_SV1, v_SV1, r_SV3, v_SV3, true);
    plot_ROE_planes_with_modes(full_times, t_orbit, d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3,SV3_modes,...
        num_orbits_modes, num_orbits_station_keep,'figures/PS5/ROE_planes_modes.png', 'figures/PS5/ROE_over_time_modes.png', 'figures/PS5/ROE_error_over_time_modes.png');
    plot_delta_v_timeline(delta_v_times, delta_v_vals, t_orbit, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
        'figures/PS5/delta_v_timeline_modes.png');

end
