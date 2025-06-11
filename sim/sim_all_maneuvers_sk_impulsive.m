function sim_all_maneuvers_sk_impulsive(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, SV1_OE_init, SV2_state_init, SV3_state_init, ...
                                     SV1_state_init, t_orbit, t_series, fig_path, title_str)
    full_times = t_series(:);
    a_chief = SV1_OE_init(1);

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

    % Maybe pick these differently later
    SV2_delta_di_max = 10.0/a_chief;
    SV3_delta_di_max = 3.0/a_chief;
    SV2_delta_de_max = 10.0/a_chief;
    SV3_delta_de_max = 3.0/a_chief;

    SV2_dphi = asin(SV2_delta_de_max/norm(SV2_roe_nom(3:4)));
    SV2_de_des = [SV2_roe_nom(3)*cos(SV2_dphi) - SV2_roe_nom(4)*sin(SV2_dphi);...
                  SV2_roe_nom(3)*sin(SV2_dphi) + SV2_roe_nom(4)*cos(SV2_dphi)];

    % Setting initial states
    SV1_state = SV1_state_init;
    SV2_state = SV2_state_init;
    SV3_state = SV3_state_init;

    state_SV3_all(1, :) = SV3_state';
    state_SV2_all(1, :) = SV2_state';
    state_SV1_all(1, :) = SV1_state';

    dt = full_times(2) - full_times(1);
    n_steps = length(full_times);

    % Initialize delta-v storage
    SV2_dv_times = [];
    SV2_dv_vals = [];
    SV2_dv_vals_complete = [];
    SV2_dv_times_complete = [];

    SV3_dv_times = [];
    SV3_dv_vals = [];
    SV3_dv_vals_complete = [];
    SV3_dv_times_complete = [];

    withJ2 = true;

    % Flag that says whether to station keep or not
    STATION_KEEP = false;
    MANEUVER_PLANNED = false;

    for i=2:n_steps
        t = full_times(i-1);

        % PRIMARY PROPAGATION
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
        [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3] = ECI2ROE_array_mean(SV1_state(1:3)', SV1_state(4:6)', SV3_state(1:3)', SV3_state(4:6)', true);
        [d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2] = ECI2ROE_array_mean(SV1_state(1:3)', SV1_state(4:6)', SV2_state(1:3)', SV2_state(4:6)', true);

        [a,e,inc,RAAN,omega,nu,M] = ECI2OE(SV1_state(1:3), SV1_state(4:6));
        SV1_oe = [a,e,inc,RAAN,omega,nu,M];
        SV1_u_val = mod(omega + M, 360);  % u = omega + M, in degrees
        SV2_roe = [d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2]/a_chief;
        SV3_roe = [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3]/a_chief;

        % This is primarily for SV3, SV2 just continuously station keeps
        if t < switch_times(1)
            % do nothing, already at initial situation
        elseif t < switch_times(2)
            MANEUVER_PLANNED = false;
            STATION_KEEP = true;
            % Calculate offset for station keeping
            SV3_roe_nom = SV3_roe_nom_mode1;
        elseif t < switch_times(3)
            if MANEUVER_PLANNED == false
                SV3_dv_vals = [];
                SV3_dv_times = [];
                [delta_v_vals, delta_v_times] = mode2_control(SV3_roe, SV3_roe_nom_mode2, t, switch_times(3), SV1_oe, SV1_u_val); % For SV3
                delta_v_idxs_sk = zeros(size(delta_v_times));
                for j = 1:length(delta_v_times)
                    [~, delta_v_idxs_sk(j)] = min(abs(full_times - delta_v_times(j)));
                end
                SV3_dv_vals = [SV3_dv_vals; delta_v_vals];
                SV3_dv_times = [SV3_dv_times, delta_v_idxs_sk];
                MANEUVER_PLANNED = true;
                STATION_KEEP = false;
            end
        elseif t < switch_times(4)
            MANEUVER_PLANNED = false;
            STATION_KEEP = true;
            % Calculate offset for station keeping
            SV3_roe_nom = SV3_roe_nom_mode2;
        elseif t < switch_times(5)
            if MANEUVER_PLANNED == false
                SV3_dv_vals = [];
                SV3_dv_times = [];
                [delta_v_vals, delta_v_times] = mode3_control(SV3_roe, SV3_roe_nom_mode3, t,switch_times(5), SV1_oe, SV1_u_val);
                delta_v_idxs_sk = zeros(size(delta_v_times));
                for j = 1:length(delta_v_times)
                    [~, delta_v_idxs_sk(j)] = min(abs(full_times - delta_v_times(j)));
                end
                SV3_dv_vals = [SV3_dv_vals; delta_v_vals];
                SV3_dv_times = [SV3_dv_times, delta_v_idxs_sk];
                MANEUVER_PLANNED = true;
                STATION_KEEP = false;
            end
        elseif t < switch_times(6)
            MANEUVER_PLANNED = false;
            STATION_KEEP = true;
            SV3_roe_nom = SV3_roe_nom_mode3;
        elseif t < switch_times(7)
            if MANEUVER_PLANNED == false
                SV3_dv_vals = [];
                SV3_dv_times = [];
                [delta_v_vals, delta_v_times] = mode4_control(SV3_roe, SV3_roe_nom_mode4, t, switch_times(7), SV1_oe, SV1_u_val);
                delta_v_idxs_sk = zeros(size(delta_v_times));
                for j = 1:length(delta_v_times)
                    [~, delta_v_idxs_sk(j)] = min(abs(full_times - delta_v_times(j)));
                end
                SV3_dv_vals = [SV3_dv_vals; delta_v_vals];
                SV3_dv_times = [SV3_dv_times, delta_v_idxs_sk];
                MANEUVER_PLANNED = true;
                STATION_KEEP = false;
            end
        elseif t < switch_times(8)
            MANEUVER_PLANNED = false;
            STATION_KEEP = true;
            % Calculate offset for station keeping
            SV3_roe_nom = SV3_roe_nom_mode4;
        else
            % sim complete
        end

        % Calculate offsets for station keeping SV3
        SV3_dphi = asin(SV3_delta_de_max/norm(SV3_roe_nom(3:4)));
        SV3_de_des = [SV3_roe_nom(3)*cos(SV3_dphi) - SV3_roe_nom(4)*sin(SV3_dphi);...
              SV3_roe_nom(3)*sin(SV3_dphi) + SV3_roe_nom(4)*cos(SV3_dphi)];

        % STATION KEEPING SV2 always
        if isempty(SV2_dv_vals)
            % SV2 eccentricity station keeping
            if norm(SV2_roe(3:4) - SV2_roe_nom(3:4)) > SV2_delta_de_max
                SV2_des_roe = SV2_roe_nom;
                SV2_des_roe(3:4) = SV2_de_des;
                % perform station keeping
                t_vals = linspace(t, t+0.5*t_orbit, 4);
                t_vals = t_vals(2:end);
                delta_v_vals_sk = naive_least_squares(t_vals, SV2_roe, SV2_des_roe, SV1_oe, SV1_u_val, t, t+0.5*t_orbit);
                delta_v_idxs_sk = zeros(size(t_vals));
                for j = 1:length(t_vals)
                    [~, delta_v_idxs_sk(j)] = min(abs(full_times - t_vals(j)));
                end
                SV2_dv_vals = [SV2_dv_vals; delta_v_vals_sk];
                SV2_dv_times = [SV2_dv_times, delta_v_idxs_sk];
            end
    
            % SV2 inclination station keeping
            if norm(SV2_roe(5:6) - SV2_roe_nom(5:6)) > SV2_delta_di_max
                SV2_des_roe = SV2_roe_nom;
                SV2_di_des = [SV2_roe_nom(5);SV2_roe_nom(6) - sign(SV2_roe(5))*SV2_delta_di_max];
                SV2_des_roe(5:6) = SV2_di_des;
                % perform station keeping 
                t_vals = linspace(t, t+0.5*t_orbit, 2);
                t_vals = t_vals(2:end);
                delta_v_vals_sk = naive_least_squares(t_vals, SV2_roe, SV2_des_roe, SV1_oe, SV1_u_val, t, t+0.5*t_orbit)';
                delta_v_idxs_sk = zeros(size(t_vals));
                for j = 1:length(t_vals)
                    [~, delta_v_idxs_sk(j)] = min(abs(full_times - t_vals(j)));
                end
                SV2_dv_vals = [SV2_dv_vals; delta_v_vals_sk];
                SV2_dv_times = [SV2_dv_times, delta_v_idxs_sk];
            end
        end
    
         % SV3 only station keeps when required
        if STATION_KEEP == true
            % SV3 eccentricity station keeping
            if isempty(SV3_dv_vals)
                if norm(SV3_roe(3:4) - SV3_roe_nom(3:4)) > SV3_delta_de_max
                    SV3_des_roe = SV3_roe_nom;
                    SV3_des_roe(3:4) = SV3_de_des;
                    % perform station keeping 
                    t_vals = linspace(t, t+0.5*t_orbit, 4);
                    t_vals = t_vals(2:end);
                    delta_v_vals_sk = naive_least_squares(t_vals, SV3_roe, SV3_des_roe, SV1_oe, SV1_u_val, t, t+0.5*t_orbit);
                    delta_v_idxs_sk = zeros(size(t_vals));
                    for j = 1:length(t_vals)
                        [~, delta_v_idxs_sk(j)] = min(abs(full_times - t_vals(j)));
                    end
                    SV3_dv_vals = [SV3_dv_vals; delta_v_vals_sk];
                    SV3_dv_times = [SV3_dv_times, delta_v_idxs_sk];
                end
        
                % SV3 inclination station keeping
                if norm(SV3_roe(5:6) - SV3_roe_nom(5:6)) > SV3_delta_di_max
                    SV3_des_roe = SV3_roe_nom;
                    SV3_di_des = [SV3_roe_nom(5);SV3_roe_nom(6) - sign(SV3_roe(5))*SV3_delta_di_max];
                    SV3_des_roe(5:6) = SV3_di_des;
                    % perform station keeping 
                    t_vals = linspace(t, t+0.5*t_orbit, 2);
                    t_vals = t_vals(2:end);
                    delta_v_vals_sk = naive_least_squares(t_vals, SV3_roe, SV3_des_roe, SV1_oe, SV1_u_val, t, t+0.5*t_orbit)';
                    delta_v_idxs_sk = zeros(size(t_vals));
                    for j = 1:length(t_vals)
                        [~, delta_v_idxs_sk(j)] = min(abs(full_times - t_vals(j)));
                    end
                    SV3_dv_vals = [SV3_dv_vals; delta_v_vals_sk];
                    SV3_dv_times = [SV3_dv_times, delta_v_idxs_sk];
                end
            end
        end

        %%% APPLY DELTA V
        SV2_time_match = find(SV2_dv_times == i, 1);
        if ~isempty(SV2_time_match)
            SV2_state(4:6) = SV2_state(4:6) + dv_RTN2ECI(SV1_state(1:3), SV1_state(4:6), SV2_dv_vals(SV2_time_match, :)'/1e3);
            SV2_dv_vals_complete = [SV2_dv_vals_complete; SV2_dv_vals(SV2_time_match, :)/1e3];
            SV2_dv_times_complete = [SV2_dv_times_complete, SV2_dv_times(SV2_time_match)];
            SV2_dv_vals(SV2_time_match, :) = [];
            SV2_dv_times(SV2_time_match) = [];
        end
        
        SV3_time_match = find(SV3_dv_times == i);
        if ~isempty(SV3_time_match)
            for idx = fliplr(SV3_time_match)  % Use fliplr to avoid indexing issues when deleting
                SV3_state(4:6) = SV3_state(4:6) + dv_RTN2ECI(SV1_state(1:3), SV1_state(4:6), SV3_dv_vals(idx, :)'/1e3);
                SV3_dv_vals_complete = [SV3_dv_vals_complete; SV3_dv_vals(idx, :)/1e3];
                SV3_dv_times_complete = [SV3_dv_times_complete, SV3_dv_times(idx)];
                SV3_dv_vals(idx, :) = [];
                SV3_dv_times(idx) = [];
            end
        end
        % if ~isempty(SV3_time_match)
        %     SV3_state(4:6) = SV3_state(4:6) + dv_RTN2ECI(SV1_state(1:3), SV1_state(4:6), SV3_dv_vals(SV3_time_match, :)'/1e3);
        %     SV3_dv_vals_complete = [SV3_dv_vals_complete; SV3_dv_vals(SV3_time_match, :)];
        %     SV3_dv_times_complete = [SV3_dv_times_complete, SV3_dv_times(SV3_time_match)];
        %     SV3_dv_vals(SV3_time_match, :) = [];
        %     SV3_dv_times(SV3_time_match) = [];
        % end

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
    plot_3D_rel_orbit(SV2_rel_pos, SV3_rel_pos, true, 'All maneuvers 3D path', 'figures/PS5/ROE_3d_all_maneuvers.png');
    % create_roe_3d_animation(SV2_rel_pos, SV3_rel_pos, 'SV2 and SV3 3D Relative Orbits (Animated)');
    
    RT = true; RN = true; NT = true;
    plot_RT_RN_projections_separate(SV2_rel_pos, SV3_rel_pos, RT, RN, NT, title_str, fig_path);
    [d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2] = ECI2ROE_array_mean(r_SV1, v_SV1, r_SV2, v_SV2, true);
    [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3] = ECI2ROE_array_mean(r_SV1, v_SV1, r_SV3, v_SV3, true);
    plot_ROE_planes(full_times, t_orbit, d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3, 'figures/PS5/ROE_planes_SV3.png', 'figures/PS5/ROE_over_time_SV3.png');
    plot_ROE_planes(full_times, t_orbit, d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2, 'figures/PS5/ROE_planes_SV2.png', 'figures/PS5/ROE_over_time_SV2.png');

    plot_ROE_planes_with_modes(full_times, t_orbit, d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2,SV2_modes,...
        num_orbits_modes, num_orbits_station_keep,'figures/PS5/ROE_planes_modes_SV2.png', 'figures/PS5/ROE_over_time_modes_SV2.png', 'figures/PS5/ROE_error_over_time_modes_SV2.png');
    plot_ROE_planes_with_modes(full_times, t_orbit, d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3,SV3_modes,...
        num_orbits_modes, num_orbits_station_keep,'figures/PS5/ROE_planes_modes_SV3.png', 'figures/PS5/ROE_over_time_modes_SV3.png', 'figures/PS5/ROE_error_over_time_modes_SV3.png');
    plot_delta_v_timeline(full_times(SV2_dv_times_complete), 1e3*SV2_dv_vals_complete, t_orbit, SV2_modes, num_orbits_modes, num_orbits_station_keep, ...
        false, 'figures/PS5/delta_v_timeline_modes_SV2.png');
    plot_delta_v_timeline(full_times(SV3_dv_times_complete), 1e3*SV3_dv_vals_complete, t_orbit, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
        false, 'figures/PS5/delta_v_timeline_modes_SV3.png');
    plot_delta_v_cumulative_timeline(full_times(SV3_dv_times_complete), SV3_dv_vals_complete*1e3, t_orbit, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
       true, 'figures/PS5/delta_v_cumulative_timeline_modes_SV3.png')
end