function sim_and_plot_control_modes_old(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
                                     SV1_OE_init, u_SV1_init, SV2_init_state, SV3_init_state, ...
                                     SV1_init_state, t_orbit, t_series, tint)
    % Constants
    num_modes = size(SV3_modes, 1);
    a_chief = SV1_OE_init(1);

    % Initialize
    SV2_roe_curr = SV2_modes(1, :)/a_chief; % divde by semi-major axis to unscale
    SV3_roe_curr = SV3_modes(1, :)/a_chief;
    t_start = 0;

    % Initialize delta-v storage
    SV2_dv_times = [];
    SV2_dv_vals = [];

    SV3_dv_times = [];
    SV3_dv_vals = [];
    dv_lower_bounds = [];

    % Loop through modes
    for mode = 2:4 % 1:num_modes
        % Determine final ROE for this mode
        SV2_roe_next = SV2_modes(mode, :)/a_chief;
        SV3_roe_next = SV3_modes(mode, :)/a_chief;

        if t_start ~= 0
            % Propagate chief from t = 0 to t_start
            t_series_mode = 0:tint:t_start;
            [~, chief_segment] = rk4_eom_ECI(t_series_mode, SV1_init_state, true);
            r_chief = chief_segment(end, 1:3)';
            v_chief = chief_segment(end, 4:6)';
            
            % Convert to orbital elements
            [a, e, i, RAAN, w, nu, M] = ECI2OE(r_chief, v_chief);
            SV1_OE_mode = [a, e, i, RAAN, w, nu];
            u_SV1_mode = mod(w + M, 360);  % u = omega + M, in degrees
        else
            SV1_OE_mode = SV1_OE_init;
            u_SV1_mode = u_SV1_init;  % u = omega + M, in degrees
        end

        % Compute mode time window
        t_mode_duration = num_orbits_modes(mode) * t_orbit;
        t_end = t_start + t_mode_duration;

        % Dynamically call modeX_control
        mode_func = str2func(sprintf('mode%d_control', mode)); % change to mode once we implement
        %[dv_SV2, times_SV2] = mode_func(SV2_roe_curr, SV2_roe_next, t_start, t_end, SV1_OE_init, u_SV1_init);
        [dv_SV3, times_SV3] = mode2_control(SV3_roe_curr, SV3_roe_next, t_start, t_end, SV1_OE_mode, u_SV1_mode);

        % Compute delta-ROEs
        d_dROE = SV3_roe_next - SV3_roe_curr;
        d_delta_e = d_dROE(3:4);
        d_delta_i = d_dROE(5:6);

        % Lower bound estimate
        dv_lb = compute_dv_lower_bound(a_chief, d_delta_e, d_delta_i);
        dv_lower_bounds = [dv_lower_bounds; dv_lb];

        fprintf('Mode %d: Total dv = %.4f m/s, Lower Bound = %.4f m/s\n', mode, sum(abs(dv_SV3), 'all'), dv_lb);
        
        % Append
        %SV2_dv_vals = [SV2_dv_vals; dv_SV2];
        %SV2_dv_times = [SV2_dv_times; times_SV2];

        SV3_dv_vals = [SV3_dv_vals; dv_SV3];
        SV3_dv_times = [SV3_dv_times, times_SV3];

        % Update current ROE
        %SV2_roe_curr = SV2_roe_next;
        SV3_roe_curr = SV3_roe_next;
        t_start = t_end;

        % Station-keeping (except after last mode)
        if mode < num_modes
            t_sk_duration = num_orbits_station_keep(mode) * t_orbit;
            t_end_sk = t_start + t_sk_duration;

            % Call station-keeping control
            %[dv_SV2_sk, times_SV2_sk] = station_keeping_control(SV2_roe_curr, SV2_roe_curr, t_start, t_end_sk, SV1_OE_init, u_SV1_init);
            %[dv_SV3_sk, times_SV3_sk] = station_keeping_control(SV3_roe_curr, SV3_roe_curr, t_start, t_end_sk, SV1_OE_init, u_SV1_init);

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
        SV3_dv_times, SV3_dv_vals, SV3_modes, num_orbits_modes, num_orbits_station_keep);  % Only SV3 gets delta-vs here, modify if SV2 gets maneuvered

    

end