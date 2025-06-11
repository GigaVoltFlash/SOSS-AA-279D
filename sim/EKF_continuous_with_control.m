function EKF_continuous_with_control(SV1_OE_init, SV2_ROE_init, SV3_ROE_init, SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, num_orbits_total, t_orbit, t_series)
    % Inputs: mean OE, mean ROEs, times

    a_chief = SV1_OE_init(1);

    % Factors for Lyaponuv feedback gain matrix  
    N = 14;
    k = 1000;

    % Settings for more things to test
    DYNAMIC_Q = false;
    ACTUATOR_NOISE = false;

    % Initialize times
    full_times = t_series(:);
    dt = full_times(2) - full_times(1);
    n_steps = length(full_times);
    switch_times = t_orbit*cumsum(reshape([num_orbits_modes; num_orbits_station_keep], 8, 1));

    % Nominal ROE for particular modes
    
    SV2_roe_nom = SV2_modes(1, :)/a_chief;
    SV3_roe_nom_mode1 = SV3_modes(1, :)/a_chief;
    SV3_roe_nom_mode2 = SV3_modes(2, :)/a_chief;
    SV3_roe_nom_mode3 = SV3_modes(3, :)/a_chief;
    SV3_roe_nom_mode4 = SV3_modes(4, :)/a_chief;
    SV3_roe_nom = SV3_roe_nom_mode1;

    % Initialize with same size as full time vector
    state_SV1_all = zeros(length(t_series), 6);
    state_SV2_all = zeros(length(t_series), 6);
    state_SV3_all = zeros(length(t_series), 6);
    OE_state_SV1_all = zeros(length(t_series), 6);
    OE_state_SV2_all = zeros(length(t_series), 6);
    OE_state_SV3_all = zeros(length(t_series), 6);

    % Initialize delta-v and acceleration storage
    SV2_dv_vals = zeros(length(t_series), 3);
    SV3_dv_vals = zeros(length(t_series), 3);
    SV2_a_vals = zeros(length(t_series), 3);
    SV3_a_vals = zeros(length(t_series), 3);
    
    % EKF Values
    x_EKF_SV2_all = zeros(length(t_series), 6);
    P_EKF_SV2_all = zeros(length(t_series), 6, 6);
    y_pred_EKF_SV2_all = zeros(length(t_series),6);
    y_post_EKF_SV2_all = zeros(length(t_series),6);
    y_actual_EKF_SV2_all = zeros(length(t_series), 6);

    x_EKF_SV3_all = zeros(length(t_series), 6);
    P_EKF_SV3_all = zeros(length(t_series), 6, 6);
    y_pred_EKF_SV3_all = zeros(length(t_series),6);
    y_post_EKF_SV3_all = zeros(length(t_series),6);
    y_actual_EKF_SV3_all = zeros(length(t_series), 6);
    
    noise_SV2_all = zeros(length(t_series), 6);
    noise_SV3_all = zeros(length(t_series), 6);

    % Use ROEs to get initial OEs for SV2 and SV3
    % Then use these OEs to get ECI positions

    a_1_init = SV1_OE_init(1);
    e_1_init = SV1_OE_init(2);
    i_1_init = SV1_OE_init(3);
    RAAN_1_init = SV1_OE_init(4);
    w_1_init = SV1_OE_init(5);
    M_1_init = SV1_OE_init(6);
    nu_1_init = rad2deg(mean2true(deg2rad(M_1_init),e_1_init));

    d_a_SV2_init = SV2_ROE_init(1);
    d_lambda_SV2_init = SV2_ROE_init(2);
    d_e_x_SV2_init = SV2_ROE_init(3);
    d_e_y_SV2_init = SV2_ROE_init(4);
    d_i_x_SV2_init = SV2_ROE_init(5);
    d_i_y_SV2_init = SV2_ROE_init(6);

    d_a_SV3_init = SV3_ROE_init(1);
    d_lambda_SV3_init = SV3_ROE_init(2);
    d_e_x_SV3_init = SV3_ROE_init(3);
    d_e_y_SV3_init = SV3_ROE_init(4);
    d_i_x_SV3_init = SV3_ROE_init(5);
    d_i_y_SV3_init = SV3_ROE_init(6);

    [a_2_init, e_2_init, i_2_init, RAAN_2_init, w_2_init, nu_2_init] = ROE2OE(a_1_init,e_1_init,i_1_init,RAAN_1_init,w_1_init,nu_1_init, ...
    d_a_SV2_init,d_lambda_SV2_init,d_e_x_SV2_init,d_e_y_SV2_init,d_i_x_SV2_init,d_i_y_SV2_init);
    [a_3_init, e_3_init, i_3_init, RAAN_3_init, w_3_init, nu_3_init] = ROE2OE(a_1_init,e_1_init,i_1_init,RAAN_1_init,w_1_init,nu_1_init, ...
    d_a_SV3_init,d_lambda_SV3_init,d_e_x_SV3_init,d_e_y_SV3_init,d_i_x_SV3_init,d_i_y_SV3_init);

    [r_ECI_SV1_init,v_ECI_SV1_init] = OE2ECI(a_1_init,e_1_init,i_1_init,RAAN_1_init,w_1_init,nu_1_init);
    [r_ECI_SV2_init,v_ECI_SV2_init] = OE2ECI(a_2_init,e_2_init,i_2_init,RAAN_2_init,w_2_init,nu_2_init);
    [r_ECI_SV3_init,v_ECI_SV3_init] = OE2ECI(a_3_init,e_3_init,i_3_init,RAAN_3_init,w_3_init,nu_3_init);

    [a_qns_SV1, e_x_SV1, e_y_SV1, i_qns_SV1, RAAN_qns_SV1, u_SV1] = ...
    OE2quasi_nonsing(a_1_init,e_1_init,i_1_init,RAAN_1_init,w_1_init,M_1_init);
    SV1_OE_state = [a_qns_SV1, e_x_SV1, e_y_SV1, i_qns_SV1, RAAN_qns_SV1, u_SV1];

    M_2_init = rad2deg(true2mean(deg2rad(nu_2_init),e_2_init));
    [a_qns_SV2, e_x_SV2, e_y_SV2, i_qns_SV2, RAAN_qns_SV2, u_SV2] = ...
    OE2quasi_nonsing(a_2_init,e_2_init,i_2_init,RAAN_2_init,w_2_init,M_2_init);
    SV2_OE_state = [a_qns_SV2, e_x_SV2, e_y_SV2, i_qns_SV2, RAAN_qns_SV2, u_SV2]; % qns

    M_3_init = rad2deg(true2mean(deg2rad(nu_3_init),e_3_init));
    [a_qns_SV3, e_x_SV3, e_y_SV3, i_qns_SV3, RAAN_qns_SV3, u_SV3] = ...
    OE2quasi_nonsing(a_3_init,e_3_init,i_3_init,RAAN_3_init,w_3_init,M_3_init);
    SV3_OE_state = [a_qns_SV3, e_x_SV3, e_y_SV3, i_qns_SV3, RAAN_qns_SV3, u_SV3]; % qns
    
    SV1_state = [r_ECI_SV1_init', v_ECI_SV1_init']';
    SV2_state = [r_ECI_SV2_init', v_ECI_SV2_init']';
    SV3_state = [r_ECI_SV3_init', v_ECI_SV3_init']';

    state_SV3_all(1, :) = SV3_state';
    state_SV2_all(1, :) = SV2_state';
    state_SV1_all(1, :) = SV1_state';
    OE_state_SV3_all(1, :) = SV3_OE_state';
    OE_state_SV2_all(1, :) = SV2_OE_state';
    OE_state_SV1_all(1, :) = SV1_OE_state';
    

    % Initial estimate and covariance for EKF (assuming perfect SV1
    % knowledge for now)
    estimate_sigma = (0.1)^2 *eye(6);
    estimate_noise = sqrtm(estimate_sigma)*randn(6,1);
    x_EKF_SV3_all(1,:) = SV3_ROE_init + 10*estimate_noise';
    P_EKF_SV3_all(1,:,:) = 100*estimate_sigma;
    x_EKF_SV2_all(1,:) = SV2_ROE_init + 10*estimate_noise';
    P_EKF_SV2_all(1,:,:) = 100*estimate_sigma;
    
    x_update_EKF_SV2 = x_EKF_SV2_all(1,:)';
    P_update_EKF_SV2 = squeeze(P_EKF_SV2_all(1,:,:));
    x_update_EKF_SV3 = x_EKF_SV3_all(1,:)';
    P_update_EKF_SV3 = squeeze(P_EKF_SV3_all(1,:,:));

    % Initial measurements and definition of noise for EKF
    RTN_sigma = (0.001)^2 * eye(3); % 0.1 m noise in each direction (defined in km)
    ECI_sigma = (0.1)^2 * eye(3); % 10 m noise in each direction (defined in km)

    [rho2, ~] = ECI2RTN_rel(SV1_state(1:3)', SV1_state(4:6)', SV2_state(1:3)', SV2_state(4:6)');
    SV2_RTN_pos = rho2';
    SV2_ECI_pos = SV2_state(1:3);

    [rho3, ~] = ECI2RTN_rel(SV1_state(1:3)', SV1_state(4:6)', SV3_state(1:3)', SV3_state(4:6)');
    SV3_RTN_pos = rho3';
    SV3_ECI_pos = SV3_state(1:3);
    
    RTN_noise_SV2 = sqrtm(RTN_sigma)*randn(3,1);
    ECI_noise_SV2 = sqrtm(ECI_sigma)*randn(3,1);
    RTN_noise_SV3 = sqrtm(RTN_sigma)*randn(3,1);
    ECI_noise_SV3 = sqrtm(ECI_sigma)*randn(3,1);

    y_actual_EKF_SV2_all(1,1:3) = SV2_RTN_pos + RTN_noise_SV2;
    y_actual_EKF_SV2_all(1,4:6) = SV2_ECI_pos + ECI_noise_SV2;
    y_actual_EKF_SV3_all(1,1:3) = SV3_RTN_pos + RTN_noise_SV3;
    y_actual_EKF_SV3_all(1,4:6) = SV3_ECI_pos + ECI_noise_SV3;

    y_pred_EKF_SV2_all(1,:) = y_actual_EKF_SV2_all(1,:);
    y_pred_EKF_SV3_all(1,:) = y_actual_EKF_SV3_all(1,:);
    y_post_EKF_SV2_all(1,:) = y_actual_EKF_SV2_all(1,:);
    y_post_EKF_SV3_all(1,:) = y_actual_EKF_SV3_all(1,:);

    % Process and measurement noise covariances
    Q3_1 = diag([1, 1, 1, 1, 1, 1]); % 0.01, 0.01
    Q3_2 = 100*diag([1, 1, 1, 1, 1, 1]); % 0.01, 0.01
    Q2 = 1e-3*diag([1, 1, 1, 1, 1, 1]); % 0.01, 0.01
    Q3 = Q3_1;

    R = 1.0*blkdiag(RTN_sigma, ECI_sigma);  % diagonal matrix with elements equal to the varianace of each measurement

    IMPULSE_MANEUVER_PLANNED = false;

    for i=2:n_steps
        t = full_times(i-1);

        %%% Planning
        % Get absolute singular chief elements from current qns
        [a_SV1,e_SV1,inc_SV1,RAAN_SV1,w_SV1,nu_SV1, M_deg_SV1] = quasi_nonsing2OE(SV1_OE_state(1), SV1_OE_state(2), SV1_OE_state(3), SV1_OE_state(4), SV1_OE_state(5), SV1_OE_state(6));
        SV1_oe = [a_SV1,e_SV1,inc_SV1,RAAN_SV1,w_SV1,nu_SV1,M_deg_SV1];
        u_deg_SV1 = mod(w_SV1 + M_deg_SV1, 360);  % u = omega + M, in degrees

        % SV2 control inputs are always determined the same way, with
        % continuous lyapunov control
        SV2_a_vals(i,:) = Lyapunov_feedback_control(x_EKF_SV2_all(i-1, :)/a_chief, SV2_roe_nom, SV1_oe, N, k, num_orbits_total);
        SV2_dv_vals(i,:) = actuator_model(dt*SV2_a_vals(i, :)/1e3, ACTUATOR_NOISE); % m/s^2 --> km/s 

        % This is primarily for SV3, SV2 just continuously station keeps     
        if t < switch_times(1)
            % do nothing, already at initial situation
        elseif t < switch_times(3) && t > switch_times(2)
            SV3_roe_nom = SV3_roe_nom_mode2; % Switches it for maneuver and station keeping
            if IMPULSE_MANEUVER_PLANNED == false
                [delta_v_vals, delta_v_times] = mode2_control(x_EKF_SV3_all(i-1, :)/a_chief, SV3_roe_nom_mode2, t, switch_times(3), SV1_oe, u_deg_SV1); % For SV3
                for j = 1:length(delta_v_times)
                    [~, idx] = min(abs(full_times - delta_v_times(j)));
                    SV3_dv_vals(idx, :) = SV3_dv_vals(idx, :) + actuator_model(delta_v_vals(j, :)/1e3, ACTUATOR_NOISE);  % m/s --> km/s
                end
                IMPULSE_MANEUVER_PLANNED = true;
            end
        elseif t < switch_times(5) && t > switch_times(4)
            SV3_roe_nom = SV3_roe_nom_mode3; % Switches it for maneuver and station keeping
            if IMPULSE_MANEUVER_PLANNED == false
                [delta_v_vals, delta_v_times] = mode3_control(x_EKF_SV3_all(i-1, :)/a_chief, SV3_roe_nom_mode3, t, switch_times(5), SV1_oe, u_deg_SV1); % For SV3
                for j = 1:length(delta_v_times)
                    [~, idx] = min(abs(full_times - delta_v_times(j)));
                    SV3_dv_vals(idx, :) = SV3_dv_vals(idx, :) + actuator_model(delta_v_vals(j, :)/1e3, ACTUATOR_NOISE);  % m/s --> km/s
                end
                IMPULSE_MANEUVER_PLANNED = true;
            end
        elseif t < switch_times(7) && t > switch_times(6)
            SV3_roe_nom = SV3_roe_nom_mode4; % Switches it for maneuver and station keeping
            if IMPULSE_MANEUVER_PLANNED == false
                [delta_v_vals, delta_v_times] = mode4_control(x_EKF_SV3_all(i-1, :)/a_chief, SV3_roe_nom_mode4, t, switch_times(7), SV1_oe, u_deg_SV1); % For SV3
                for j = 1:length(delta_v_times)
                    [~, idx] = min(abs(full_times - delta_v_times(j)));
                    SV3_dv_vals(idx, :) = SV3_dv_vals(idx, :) + actuator_model(delta_v_vals(j, :)/1e3, ACTUATOR_NOISE);  % m/s --> km/s
                end
                IMPULSE_MANEUVER_PLANNED = true;
            end        
        else
            % Run station keeping 
            IMPULSE_MANEUVER_PLANNED = false;
            SV3_a_vals(i, :) = Lyapunov_feedback_control(x_EKF_SV3_all(i-1, :)/a_chief, SV3_roe_nom, SV1_oe, N, k, num_orbits_station_keep(2)); % Currently hardcoded number of orbits for station keep
            SV3_dv_vals(i, :) = actuator_model(dt*SV3_a_vals(i, :)/1e3, ACTUATOR_NOISE); % m/s^2 --> m/s
        end

        %%% PROPAGATE GROUND TRUTH
        % Use GVE propagated for ground truth qns mean OE, then convert to
        % ground truth ECI positions
        SV1_OE_state_inter = SV1_OE_state + (dt*secular_J2_and_drag(t, SV1_OE_state))'; % propagate chief qns OE using GVE
        SV1_OE_state = wrap_QNSOE(SV1_OE_state_inter); % wraps u to be within 0-360
        SV2_OE_state_inter = SV2_OE_state + (dt*secular_J2_and_drag(t, SV2_OE_state))'; % propagate SV2 qns OE using GVE (take out later)
        SV2_OE_state_pre_ctrl = wrap_QNSOE(SV2_OE_state_inter);
        SV3_OE_state_inter = SV3_OE_state + (dt*secular_J2_and_drag(t, SV3_OE_state))'; % propagate SV3 qns OE using GVE (take out later)
        SV3_OE_state_pre_ctrl = wrap_QNSOE(SV3_OE_state_inter);

        %%% APPLY CONTROL
        % Convert mean OE to ECI
        [SV1_state_pre_ctrl, SV2_state_pre_ctrl, SV3_state_pre_ctrl, ...
          r_ECI_SV1_pre_ctrl, v_ECI_SV1_pre_ctrl, r_ECI_SV2_pre_ctrl, v_ECI_SV2_pre_ctrl, r_ECI_SV3_pre_ctrl, v_ECI_SV3_pre_ctrl] = ...
          chief_deputy_OEs_qns2ECIs(SV1_OE_state, SV2_OE_state_pre_ctrl, SV3_OE_state_pre_ctrl);
        
        SV2_state_post_ctrl = SV2_state_pre_ctrl;
        SV3_state_post_ctrl = SV3_state_pre_ctrl;

        SV2_state_post_ctrl(4:6) = SV2_state_pre_ctrl(4:6) + dv_RTN2ECI(SV1_state(1:3), SV1_state(4:6), SV2_dv_vals(i,:)')';
        SV3_state_post_ctrl(4:6) = SV3_state_pre_ctrl(4:6) + dv_RTN2ECI(SV1_state(1:3), SV1_state(4:6), SV3_dv_vals(i,:)')';

        % Convert back to osculating OE from post-control ECI
        SV2_OE_state_osc = ECI2quasi_nonsing(SV2_state_post_ctrl);
        SV3_OE_state_osc = ECI2quasi_nonsing(SV3_state_post_ctrl);

        % Convert osculating OE back to mean OE for propagation
        SV2_OE_state = qns_osc2qns_mean(SV2_OE_state_osc);
        SV3_OE_state = qns_osc2qns_mean(SV3_OE_state_osc);
     
        % Get final ECI state of chief and deputy after control is applied
        [SV1_state, SV2_state, SV3_state, ...
          r_ECI_SV1, v_ECI_SV1, r_ECI_SV2, v_ECI_SV2, r_ECI_SV3, v_ECI_SV3] = ...
          chief_deputy_OEs_qns2ECIs(SV1_OE_state, SV2_OE_state, SV3_OE_state);

        %%% GENERATE MEASUREMENTS
        % EKF Generate Measurements
        [rho2, ~] = ECI2RTN_rel(r_ECI_SV1', v_ECI_SV1', r_ECI_SV2', v_ECI_SV2');
        SV2_RTN_pos = rho2';
        SV2_ECI_pos = SV2_state(1:3)';
        RTN_noise = sqrtm(RTN_sigma)*randn(3,1); % use sigma value defined above
        ECI_noise = sqrtm(ECI_sigma)*randn(3,1);
        SV2_RTN_measurement = SV2_RTN_pos + RTN_noise;
        SV2_ECI_measurement = SV2_ECI_pos + ECI_noise;
        noise_vector_SV2 = [RTN_noise;ECI_noise];
        y_actual_EKF_SV2 = [SV2_RTN_measurement;SV2_ECI_measurement];

        [rho3, ~] = ECI2RTN_rel(r_ECI_SV1', v_ECI_SV1', r_ECI_SV3', v_ECI_SV3');
        SV3_RTN_pos = rho3';
        SV3_ECI_pos = SV3_state(1:3)';
        RTN_noise = sqrtm(RTN_sigma)*randn(3,1); % use sigma value defined above
        ECI_noise = sqrtm(ECI_sigma)*randn(3,1);
        SV3_RTN_measurement = SV3_RTN_pos + RTN_noise;
        SV3_ECI_measurement = SV3_ECI_pos + ECI_noise;
        noise_vector_SV3 = [RTN_noise;ECI_noise];
        y_actual_EKF_SV3 = [SV3_RTN_measurement;SV3_ECI_measurement];
        
        %%% RUN EKF
        u_SV2 = SV2_dv_vals(i, :)';
        u_SV3 = SV3_dv_vals(i, :)'; 

        if DYNAMIC_Q
            if IMPULSE_MANEUVER_PLANNED==true & any(SV3_dv_vals(i, :))==true
                Q3 = Q3_2; % Use a much larger Q at the time-steps that we have an impulse
            else
                Q3 = Q3_1; 
            end
        end

        % Run EKF for SV2
        [x_update_EKF_SV2, P_update_EKF_SV2, y_pred_EKF_SV2, y_post_EKF_SV2] = ekf_roes_with_control(x_update_EKF_SV2, y_actual_EKF_SV2, P_update_EKF_SV2, SV1_state, SV1_OE_state, Q2, R, dt, u_SV2);

        % Run EKF for SV3
        [x_update_EKF_SV3, P_update_EKF_SV3, y_pred_EKF_SV3, y_post_EKF_SV3] = ekf_roes_with_control(x_update_EKF_SV3, y_actual_EKF_SV3, P_update_EKF_SV3, SV1_state, SV1_OE_state, Q3, R, dt, u_SV3);
        
        % Save states to buffer
        state_SV3_all(i, :) = SV3_state';
        state_SV2_all(i, :) = SV2_state';
        state_SV1_all(i, :) = SV1_state';
        OE_state_SV3_all(i, :) = SV3_OE_state';
        OE_state_SV2_all(i, :) = SV2_OE_state';
        OE_state_SV1_all(i, :) = SV1_OE_state';


        % EKF buffer
        y_actual_EKF_SV2_all(i,:) = y_actual_EKF_SV2;
        y_pred_EKF_SV2_all(i, :) = y_pred_EKF_SV2;
        y_post_EKF_SV2_all(i, :) = y_post_EKF_SV2;
        x_EKF_SV2_all(i, :) = x_update_EKF_SV2;
        P_EKF_SV2_all(i,:,:) = P_update_EKF_SV2;
        noise_SV2_all(i, :) = noise_vector_SV2;

        y_actual_EKF_SV3_all(i,:) = y_actual_EKF_SV3;
        y_pred_EKF_SV3_all(i, :) = y_pred_EKF_SV3;
        y_post_EKF_SV3_all(i, :) = y_post_EKF_SV3;
        x_EKF_SV3_all(i, :) = x_update_EKF_SV3;
        P_EKF_SV3_all(i,:,:) = P_update_EKF_SV3;
        noise_SV3_all(i, :) = noise_vector_SV3;
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

    % What do we need to plot?
    
        % 1. True and Estimated States with covariance overlay
        % 2. Estimation Error
        % 3. True Statistics (mean and standard deviation) (not really a plot)
        % 4. Pre-fit residual (y_act - y_pred) and post-fit residual (y_act - y_with_x_update) against input noise
    
    [d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2] = ECI2ROE_array_mean(r_SV1, v_SV1, r_SV2, v_SV2, true); 
    [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3] = ECI2ROE_array_mean(r_SV1, v_SV1, r_SV3, v_SV3, true); 

    ROE_SV3_true = [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3];
    ROE_SV2_true = [d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2];

    pre_fit_residual_SV2 = y_actual_EKF_SV2_all - y_pred_EKF_SV2_all;
    post_fit_residual_SV2 = y_actual_EKF_SV2_all - y_post_EKF_SV2_all;

    pre_fit_residual_SV3 = y_actual_EKF_SV3_all - y_pred_EKF_SV3_all;
    post_fit_residual_SV3 = y_actual_EKF_SV3_all - y_post_EKF_SV3_all;
    EKF_error_SV3 = x_EKF_SV3_all - ROE_SV3_true;
    EKF_error_SV2 = x_EKF_SV2_all - ROE_SV2_true;

    % SV2 Plots
    plot_ROE_comparison_with_cov(full_times, t_orbit, ROE_SV2_true, x_EKF_SV2_all,  P_EKF_SV2_all, 'Ground Truth', 'EKF',...
        '','figures/PS9/ROE_planes_SV2_comparison.png', '', 'figures/PS9/ROE_over_time_SV2_comparison.png');
    
    plot_EKF_error(full_times, t_orbit, EKF_error_SV2, P_EKF_SV2_all, 'figures/PS9/EKF_error_SV2.png');
    
    plot_EKF_residuals(full_times, t_orbit, pre_fit_residual_SV2, post_fit_residual_SV2, noise_SV2_all, 'figures/PS9/residuals_SV2.png');

    % SV3 Plots
    plot_ROE_comparison_with_cov(full_times, t_orbit, ROE_SV3_true, x_EKF_SV3_all,  P_EKF_SV3_all, 'Ground Truth', 'EKF',...
        '','figures/PS9/ROE_planes_SV3_comparison.png', '', 'figures/PS9/ROE_over_time_SV3_comparison.png');
    
    plot_EKF_error(full_times, t_orbit, EKF_error_SV3, P_EKF_SV3_all, 'figures/PS9/EKF_error_SV3.png');
    
    plot_EKF_residuals(full_times, t_orbit, pre_fit_residual_SV3, post_fit_residual_SV3, noise_SV3_all, 'figures/PS9/residuals_SV3.png');

    % Mode plots
    plot_ROE_planes_with_modes(full_times, t_orbit, d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2,SV2_modes,...
        num_orbits_modes, num_orbits_station_keep,'figures/PS9/ROE_planes_modes_SV2.png', 'figures/PS9/ROE_over_time_modes_SV2.png', 'figures/PS9/ROE_error_over_time_modes_SV2.png');
    plot_ROE_planes_with_modes(full_times, t_orbit, d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3,SV3_modes,...
        num_orbits_modes, num_orbits_station_keep,'figures/PS9/ROE_planes_modes_SV3.png', 'figures/PS9/ROE_over_time_modes_SV3.png', 'figures/PS9/ROE_error_over_time_modes_SV3.png');

    % Delta V plots to see control inputs
    plot_delta_v_timeline(full_times(:), SV2_dv_vals*1e3, t_orbit, SV2_modes, num_orbits_modes, num_orbits_station_keep, ...
        true, 'figures/PS9/delta_v_timeline_modes_SV2.png');
    plot_delta_v_timeline(full_times(:), SV3_dv_vals*1e3, t_orbit, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
        false, 'figures/PS9/delta_v_timeline_modes_SV3.png');
    plot_delta_v_cumulative_timeline(full_times(:), SV3_dv_vals*1e3, t_orbit, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
       true, 'figures/PS9/delta_v_cumulative_timeline_modes_SV3.png')

    % Plot the 3D orbit
    plot_3D_rel_orbit(SV2_rel_pos, SV3_rel_pos, true, 'All maneuvers 3D path', 'figures/PS9/ROE_3d_all_maneuvers.png');
    % final_EKF_mean_err_SV2 = EKF_error_SV2(end, :);
    % final_EKF_mean_err_SV3 = EKF_error_SV3(end, :);
    % final_EKF_std_err_SV2 = sqrt(diag(squeeze(P_EKF_SV2_all(end, :, :))));
    % final_EKF_std_err_SV3 = sqrt(diag(squeeze(P_EKF_SV3_all(end, :, :))));
    % disp('Final EKF Mean Error SV2:');
    % disp(final_EKF_mean_err_SV2);
    % disp('Final EKF Std Error SV2:');
    % disp(final_EKF_std_err_SV2);
    % 
    % disp('Final EKF Mean Error SV3:');
    % disp(final_EKF_mean_err_SV3);
    % disp('Final EKF Std Error SV3:');
    % disp(final_EKF_std_err_SV3);

    % Print as LaTeX matrices (3 decimal places)
    % fprintf('Final EKF Mean Error SV2 (LaTeX):\n');
    % fprintf('\\begin{bmatrix} %s \\end{bmatrix}\n', ...
    %     strjoin(arrayfun(@(x) sprintf('%.3f', x), final_EKF_mean_err_SV2, 'UniformOutput', false), ' \\\\ '));
    % fprintf('Final EKF Std Error SV2 (LaTeX):\n');
    % fprintf('\\begin{bmatrix} %s \\end{bmatrix}\n', ...
    %     strjoin(arrayfun(@(x) sprintf('%.3f', x), final_EKF_std_err_SV2, 'UniformOutput', false), ' \\\\ '));
    % 
    % fprintf('Final EKF Mean Error SV3 (LaTeX):\n');
    % fprintf('\\begin{bmatrix} %s \\end{bmatrix}\n', ...
    %     strjoin(arrayfun(@(x) sprintf('%.3f', x), final_EKF_mean_err_SV3, 'UniformOutput', false), ' \\\\ '));
    % fprintf('Final EKF Std Error SV3 (LaTeX):\n');
    % fprintf('\\begin{bmatrix} %s \\end{bmatrix}\n', ...
    %     strjoin(arrayfun(@(x) sprintf('%.3f', x), final_EKF_std_err_SV3, 'UniformOutput', false), ' \\\\ '));
    
end