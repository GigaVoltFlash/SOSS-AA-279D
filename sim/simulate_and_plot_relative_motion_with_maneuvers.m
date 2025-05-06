function simulate_and_plot_relative_motion_with_maneuvers(t, SV2_roe_init, state_abs_SV2_init, ...
    SV3_roe_init, state_abs_SV3_init, SV1_oe_init, SV1_state_init, ...
    fig_path, title_str, delta_v_times, delta_v_vals)

    % Propagate chief (with J2)
    [~, chief_state] = rk4_eom_ECI(t, SV1_state_init, true);
    r_SV1 = chief_state(:,1:3);
    v_SV1 = chief_state(:,4:6);

    % Propagate SV2 (no impulses)
    [~, state_SV2] = rk4_eom_ECI(t, state_abs_SV2_init, true);
    r_SV2 = state_SV2(:,1:3);
    v_SV2 = state_SV2(:,4:6);

    % --- Split SV3 propagation into segments with impulsive delta-v's ---
    full_times = t(:);
    man_times = [0; delta_v_times(:); full_times(end)];
    % Initialize with same size as full time vector
    state_SV3_all = zeros(length(t), 6);
    
    % Index tracker for filling in state_SV3_all
    idx_start = 1;
    state_curr = state_abs_SV3_init';
    
    % Loop over maneuver segments
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
            dv = (delta_v_vals(i-1, :)/1e3)';  % delta-v at maneuver time (m/s --> km/s)
            state_curr(4:6) = state_curr(4:6) + dv;
        end
    
        % Propagate with rk4
        [~, state_segment] = rk4_eom_ECI(t_segment, state_curr, true);
    
        % Save to full array
        state_SV3_all(idx_range, :) = state_segment;
    
        % Update for next segment
        state_curr = state_segment(end, :)';
    end

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
    plot_RT_RN_projections_separate(SV2_rel_pos, SV3_rel_pos, RT, RN, NT, title_str, fig_path);
end
