function sim_and_plot_control_modes(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
                                     SV1_OE_init, u_SV1_init, SV2_init_state, SV3_init_state, ...
                                     SV1_init_state, t_orbit, t_series, tint, fig_path, title_str)
    % Constants
    num_modes = size(SV3_modes, 1);
    a_chief = SV1_OE_init(1);

    % Initialize
    SV2_roe_curr = SV2_modes(1, :)/a_chief; % divde by semi-major axis to unscale
    SV3_roe_curr = SV3_modes(1, :)/a_chief;
    t_start = 0;
    %state_curr_SV3_all = state_abs_SV3_init';
    %state_curr_SV2_all = state_abs_SV2_init';
    %state_curr_SV1_all = SV1_state_init';

    % Initialize delta-v storage
    SV2_dv_times = [];
    SV2_dv_vals = [];

    SV3_dv_times = [];
    SV3_dv_vals = [];

    % Initialize with same size as full time vector
    state_SV1_all = zeros(length(t_series), 6);
    state_SV2_all = zeros(length(t_series), 6);
    state_SV3_all = zeros(length(t_series), 6);

    % Loop through modes
    for mode = 2:3 % 1:num_modes
        % Determine final ROE for this mode
        SV2_roe_next = SV2_modes(mode, :)/a_chief;
        SV3_roe_next = SV3_modes(mode, :)/a_chief;

        % Compute mode time window
        t_mode_duration = num_orbits_modes(mode) * t_orbit;
        t_end = t_start + t_mode_duration;
        t_series_mode = t_start:tint:t_end;

        if t_start ~= 0
            % Propagate chief from t = 0 to t_start
            %t_series_mode = 0:tint:t_start;
            %[~, chief_segment] = rk4_eom_ECI(t_series_mode, SV1_init_state, true);
            r_chief = state_curr_SV1(1:3);
            v_chief = state_curr_SV1(4:6);
            
            % Convert to orbital elements
            [a, e, i, RAAN, w, nu, M] = ECI2OE(r_chief, v_chief);
            SV1_OE_mode = [a, e, i, RAAN, w, nu];
            u_SV1_mode = mod(w + M, 360);  % u = omega + M, in degrees
        else
            SV1_OE_mode = SV1_OE_init;
            u_SV1_mode = u_SV1_init;  % u = omega + M, in degrees
        end
       
        % Dynamically call modeX_control
        %mode_func = str2func(sprintf('mode%d_control', mode)); % change to mode once we implement
        %[dv_SV2, times_SV2] = mode_func(SV2_roe_curr, SV2_roe_next, t_start, t_end, SV1_OE_init, u_SV1_init);
        [dv_SV3, times_SV3] = mode2_control(SV3_roe_curr, SV3_roe_next, t_start, t_end, SV1_OE_mode, u_SV1_mode);

        delta_v_times = times_SV3;
        delta_v_vals = dv_SV3;
        
        full_times = t_series_mode(:);
        man_times = [0; delta_v_times(:); full_times(end)];
        % Initialize with same size as man time vector
        %state_SV1_man = zeros(length(t_series_mode), 6);
        %state_SV2_man = zeros(length(t_series_mode), 6);
        %state_SV3_man = zeros(length(t_series_mode), 6);

        % Index tracker for filling in state_SV3_all
        if t_start == 0
            state_curr_SV3 = SV3_init_state';
            state_curr_SV2 = SV2_init_state';
            state_curr_SV1 = SV1_init_state';
        end

        % Loop over maneuver segments
        % Currently just maneuvering SV3
        for i = 1:length(man_times) - 1
            % Find index range for this time segment
            idx_range = find(t_series >= man_times(i) & t_series <= man_times(i+1));
        
            if isempty(idx_range)
                continue; % Skip empty segments (shouldn’t happen with valid times)
            end
        
            % Time vector for this segment
            t_segment = t_series(idx_range);
        
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
        
        % Append
        %SV2_dv_vals = [SV2_dv_vals; dv_SV2];
        %SV2_dv_times = [SV2_dv_times; times_SV2];

        %SV3_dv_vals = [SV3_dv_vals; dv_SV3];
        %SV3_dv_times = [SV3_dv_times, times_SV3];

        % Update current ROE
        %SV2_roe_curr = SV2_roe_next;
        %SV3_roe_curr = SV3_roe_next;

        r_SV1_curr = state_curr_SV1(1:3);
        v_SV1_curr = state_curr_SV1(4:6);
        r_SV3_curr = state_curr_SV3(1:3);
        v_SV3_curr = state_curr_SV3(4:6);

        [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3] = ...
            ECI2ROE_array_mean(r_SV1_curr, v_SV1_curr, r_SV3_curr, v_SV3_curr, true);
        SV3_roe_curr = [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3]/a_chief;

        t_start = t_end;

        % Station-keeping (except after last mode)
        if mode < num_modes
            t_sk_duration = num_orbits_station_keep(mode) * t_orbit;
            t_end_sk = t_start + t_sk_duration;
            %[dv_SV2_sk, times_SV2_sk] = mode_station_keeping(SV2_roe_curr, SV2_roe_curr, t_start, t_end_sk, SV1_OE_init, u_SV1_init);
            %[dv_SV3_sk, times_SV3_sk] = mode_station_keeping(SV3_roe_curr, SV3_roe_curr, t_start, t_end_sk, SV1_OE_init, u_SV1_init);

            % Append
            %SV2_dv_vals = [SV2_dv_vals; dv_SV2_sk];
            %SV2_dv_times = [SV2_dv_times; times_SV2_sk];

            %SV3_dv_vals = [SV3_dv_vals; dv_SV3_sk];
            %SV3_dv_times = [SV3_dv_times; times_SV3_sk];

            % Advance time
            t_start = t_end_sk;
        end
    end

    SV2_roe_init = SV2_modes(1,:)/a_chief; % divde by semi-major axis to unscale
    SV3_roe_init = SV3_modes(1,:)/a_chief;
    
    % Simulate and plot using full delta-v profiles
    simulate_and_plot_relative_motion_with_maneuvers(t_series, t_orbit,...
        SV2_roe_init, SV2_init_state, ...
        SV3_roe_init, SV3_init_state, ...
        SV1_OE_init, SV1_init_state, ...
        fig_path, title_str, ...
        SV3_dv_times, SV3_dv_vals, SV3_modes, num_orbits_modes);  % Only SV3 gets delta-vs here, modify if SV2 gets maneuvered

end