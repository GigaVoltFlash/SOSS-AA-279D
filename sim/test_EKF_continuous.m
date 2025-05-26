function test_EKF_continuous(SV2_modes, SV3_modes, num_orbits_modes, num_orbits_station_keep, SV1_OE_init, SV2_state_init, SV3_state_init, ...
                                     SV1_state_init, t_orbit, t_series)
    full_times = t_series(:);
    switch_times = t_orbit*cumsum(reshape([num_orbits_modes; num_orbits_station_keep], 8, 1));
    dt = full_times(2) - full_times(1);
    n_steps = length(full_times);

    a_chief = SV1_OE_init(1);

    withJ2 = true;

    % Factors for Lyaponuv feedback gain matrix  
    N = 14;
    k = 1000;

    % Initialize with same size as full time vector
    state_SV1_all = zeros(length(t_series), 6);
    state_SV2_all = zeros(length(t_series), 6);
    state_SV3_all = zeros(length(t_series), 6);

    SV1_OE_test = zeros(length(t_series), 6);

    % Initialize delta-v and acceleration storage
    SV2_dv_vals = zeros(length(t_series), 3);
    SV3_dv_vals = zeros(length(t_series), 3);
    SV2_a_vals = zeros(length(t_series), 3);
    SV3_a_vals = zeros(length(t_series), 3);
    
    
    % EKF Values
    OE_EKF_SV1_all = zeros(length(t_series), 6);
    OE_EKF_SV3_all = zeros(length(t_series), 6);
    ROE_EKF_SV3_all = zeros(length(t_series), 6);
    x_pred_EKF_SV3_all = zeros(length(t_series), 6);
    P_pred_EKF_SV3_all = zeros(length(t_series), 6, 6);
    x_update_EKF_SV3_all = zeros(length(t_series), 6);
    P_update_EKF_SV3_all = zeros(length(t_series), 6, 6);
    y_pred_EKF_SV3_all = zeros(length(t_series),6);
    y_actual_EKF_SV3_all = zeros(length(t_series), 6);

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
    
    
    % Initial states for EKF
    [a_init_SV1,e_init_SV1,i_init_SV1,RAAN_init_SV1,omega_init_SV1,nu_init_SV1] = ECI2OE(SV1_state_init(1:3), SV1_state_init(4:6));
    M_init_SV1_rad = true2mean(deg2rad(nu_init_SV1),e_init_SV1);
    M_init_SV1 = rad2deg(M_init_SV1_rad);
    [a_qns_SV1, e_x_SV1, e_y_SV1, i_qns_SV1, RAAN_qns_SV1, u_SV1] = ...
    OE2quasi_nonsing(a_init_SV1,e_init_SV1,i_init_SV1,RAAN_init_SV1,omega_init_SV1,M_init_SV1);
    SV1_OE_state = [a_qns_SV1, e_x_SV1, e_y_SV1, i_qns_SV1, RAAN_qns_SV1, u_SV1];
    OE_EKF_SV1_all(1,:) = SV1_OE_state; % qns
    
    [d_a_init_SV3, d_lambda_init_SV3, d_e_x_init_SV3, d_e_y_init_SV3, d_i_x_init_SV3, d_i_y_init_SV3] = ...
        ECI2ROE_array(SV1_state_init(1:3)', SV1_state_init(4:6)', SV3_state_init(1:3)', SV3_state_init(4:6)');
    SV3_ROE_state = [d_a_init_SV3, d_lambda_init_SV3, d_e_x_init_SV3, d_e_y_init_SV3, d_i_x_init_SV3, d_i_y_init_SV3]/a_chief; % unscaling by a_chief
    ROE_EKF_SV3_all(1, :) = SV3_ROE_state*a_chief; % storing scaled ROEs
        %SV3_ROE_prime_state
    
    [a_qns_SV3, e_x_SV3, e_y_SV3, i_qns_SV3, RAAN_qns_SV3, u_SV3] = ...
        ROE2quasi_nonsing(a_qns_SV1, e_x_SV1, e_y_SV1, i_qns_SV1, RAAN_qns_SV1, u_SV1,d_a_init_SV3, d_lambda_init_SV3, d_e_x_init_SV3, d_e_y_init_SV3, d_i_x_init_SV3, d_i_y_init_SV3);
    SV3_OE_state = [a_qns_SV3, e_x_SV3, e_y_SV3, i_qns_SV3, RAAN_qns_SV3, u_SV3];
    OE_EKF_SV3_all(1,:) = SV3_OE_state; % qns

    % Initial estimate and covariance for EKF (assuming perfect SV1
    % knowledge for now)
    estimate_sigma = 0.1*eye(6);
    %estimate_sigma(2,2) = 1e2; 
    estimate_noise = sqrtm(estimate_sigma)*randn(6,1);
    x_pred_EKF_SV3_all(1,:) = ROE_EKF_SV3_all(1, :); % + estimate_noise';
    P_pred_EKF_SV3_all(1,:,:) = estimate_sigma;
    x_update_EKF_SV3_all(1,:) = x_pred_EKF_SV3_all(1,:);
    P_update_EKF_SV3_all(1,:,:) = P_pred_EKF_SV3_all(1,:,:);

    %x_update_EKF_SV3_prime_unscaled = ROE2ROE_prime((x_update_EKF_SV3_all(1,:)/a_chief)',SV1_OE_state'); % unscaled, qns prime ROE
    %P_pred_EKF_SV3 = squeeze(P_pred_EKF_SV3_all(1,:,:));

    x_update_EKF_SV3 = x_update_EKF_SV3_all(1,:)';
    P_update_EKF_SV3 = squeeze(P_update_EKF_SV3_all(1,:,:));

    % Generates STM and unscaled roe output at every time step
    [STMs, ROEs_SV3_STM] = roe_stm_j2(t_series, SV3_ROE_state*a_chief, SV1_OE_init);
    %ROEs_SV3_STM = ROEs_SV3_STM_unscaled*a_chief;

    % Initial measurements and definition of noise for EKF
    RTN_sigma = 0.000001 * eye(3); % 1 m noise in each direction .000001
    ECI_sigma = 1 * eye(3); % 100 m noise in each direction .001 or 1

    [rho3, ~] = ECI2RTN_rel(SV1_state(1:3)', SV1_state(4:6)', SV3_state(1:3)', SV3_state(4:6)');
    SV3_RTN_pos = rho3';
    SV3_ECI_pos = SV3_state(1:3);
    
    RTN_noise = sqrtm(RTN_sigma)*randn(3,1);
    
    ECI_noise = sqrtm(ECI_sigma)*randn(3,1);

    y_actual_EKF_SV3_all(1,1:3) = SV3_RTN_pos + RTN_noise;
    y_actual_EKF_SV3_all(1,4:6) = SV3_ECI_pos + ECI_noise;

    y_pred_EKF_SV3_all(1,:) = y_actual_EKF_SV3_all(1,:);

    % Process and measurement noise covariances
    Q = 0.1* estimate_sigma; % similar to P_0 but much smaller
    %Q(1,1) = 1e4;
    %Q(2,2) = 1e3;
    %Q(5,5) = 1e1;
    %Q(6,6) = 1e1;

    R = blkdiag(RTN_sigma, ECI_sigma); % diagonal matrix with elements equal to the varianace of each measurement

    for i=2:n_steps
        t = full_times(i-1);

        % PRIMARY PROPAGATION (ECI RK4 STEPS)
        k1_SV1 = eom_ECI(t, SV1_state, withJ2);
        k2_SV1 = eom_ECI(t + dt/2, SV1_state + dt/2 * k1_SV1, withJ2);
        k3_SV1 = eom_ECI(t + dt/2, SV1_state + dt/2 * k2_SV1, withJ2);
        k4_SV1 = eom_ECI(t + dt,   SV1_state + dt   * k3_SV1, withJ2);
        %SV1_state = SV1_state + dt/6 * (k1_SV1 + 2*k2_SV1 + 2*k3_SV1 + k4_SV1);
        
        k1_SV2 = eom_ECI(t, SV2_state, withJ2);
        k2_SV2 = eom_ECI(t + dt/2, SV2_state + dt/2 * k1_SV2, withJ2);
        k3_SV2 = eom_ECI(t + dt/2, SV2_state + dt/2 * k2_SV2, withJ2);
        k4_SV2 = eom_ECI(t + dt,   SV2_state + dt   * k3_SV2, withJ2);
        SV2_state = SV2_state + dt/6 * (k1_SV2 + 2*k2_SV2 + 2*k3_SV2 + k4_SV2);
        
        k1_SV3 = eom_ECI(t, SV3_state, withJ2);
        k2_SV3 = eom_ECI(t + dt/2, SV3_state + dt/2 * k1_SV3, withJ2);
        k3_SV3 = eom_ECI(t + dt/2, SV3_state + dt/2 * k2_SV3, withJ2);
        k4_SV3 = eom_ECI(t + dt,   SV3_state + dt   * k3_SV3, withJ2);
        %SV3_state = SV3_state + dt/6 * (k1_SV3 + 2*k2_SV3 + 2*k3_SV3 + k4_SV3);

        % Funky propagation fixes
        SV1_OE_state_inter = SV1_OE_state + (dt*secular_J2(t, SV1_OE_state))'; % propagate chief qns OE using GVE
        SV1_OE_state = wrap_QNSOE(SV1_OE_state_inter);
        SV3_OE_state_inter = SV3_OE_state + (dt*secular_J2(t, SV3_OE_state))'; % propagate SV3 qns OE using GVE (take out later)
        SV3_OE_state = wrap_QNSOE(SV3_OE_state_inter);

        d_a_SV3_STM = ROEs_SV3_STM(i,1);
        d_lambda_SV3_STM =  ROEs_SV3_STM(i,2);
        d_e_x_SV3_STM = ROEs_SV3_STM(i,3);
        d_e_y_SV3_STM = ROEs_SV3_STM(i,4);
        d_i_x_SV3_STM = ROEs_SV3_STM(i,5);
        d_i_y_SV3_STM = ROEs_SV3_STM(i,6);

        a_o = SV1_OE_state(1);
        e_x_o = SV1_OE_state(2);
        e_y_o = SV1_OE_state(3);
        i_o = SV1_OE_state(4);
        RAAN_o = SV1_OE_state(5);
        u_o = SV1_OE_state(6);

        [r_ECI_SV3,v_ECI_SV3] = ROE2ECI(a_o,e_x_o,e_y_o,i_o,RAAN_o,u_o, ...
        d_a_SV3_STM,d_lambda_SV3_STM,d_e_x_SV3_STM,d_e_y_SV3_STM,d_i_x_SV3_STM,d_i_y_SV3_STM);

        [a_o,e_o,i_o,RAAN_o,w_o,nu_o, M_o] = quasi_nonsing2OE(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o); 
        [r_ECI_SV1,v_ECI_SV1] = OE2ECI(a_o,e_o,i_o,RAAN_o,w_o,nu_o);

        SV1_OE_sing(i,:) = [a_o,e_o,i_o,RAAN_o,w_o,nu_o];

        SV3_state = [r_ECI_SV3', v_ECI_SV3']';
        SV1_state = [r_ECI_SV1', v_ECI_SV1']';

        % [a,e,i,RAAN,omega,nu,M] = ECI2OE(r_IJK,v_IJK)
        %[a_test,e_test,i_test,RAAN_test,omega_test,~,M_test] = ECI2OE(SV1_state(1:3),SV1_state(4:6));
        % [a, e_x, e_y, i, RAAN, u] = OE2quasi_nonsing(a, e, i, RAAN, w, M)
        %[a_test, e_x_test, e_y_test, i_test, RAAN_test, u_test]= OE2quasi_nonsing(a_test,e_test,i_test,RAAN_test,omega_test,M_test);
        %SV1_OE_qns_test = [a_test, e_x_test, e_y_test, i_test, RAAN_test, u_test];

        [STM_big,ROE_test] = roe_stm_j2(dt, x_update_EKF_SV3', SV1_OE_sing(i,:));
        STM_curr = squeeze(STM_big);

        % EKF Mean Prediction (Chief OE, and ROE Euler Integration

        x_update_EKF_SV3_prime = ROE2ROE_prime(x_update_EKF_SV3,SV1_OE_state);
        x_update_EKF_SV3_prime_unscaled = x_update_EKF_SV3_prime/a_chief;
        
        x_update_EKF_SV3_prime_unscaled = x_update_EKF_SV3_prime_unscaled + dt*eom_ROE_prime(t, SV3_OE_state, SV1_OE_state); % propagate ROE with Euler

        x_pred_EKF_SV3_unscaled = ROE_prime2ROE(x_update_EKF_SV3_prime_unscaled,SV1_OE_state);
        x_pred_EKF_SV3 = (x_pred_EKF_SV3_unscaled')*a_chief; % scaling by a_chief

        %x_pred_EKF_SV3 = STM_curr*x_update_EKF_SV3;

        % EKF Covariance Prediction
        P_pred_EKF_SV3 = STM_curr * P_update_EKF_SV3 * STM_curr' + Q;

        % EKF Predict Measurements
        y_pred_EKF_SV3 = measurement_model_SV3(x_pred_EKF_SV3,SV1_OE_state);

        % EKF Generate Measurements
        [rho3, ~] = ECI2RTN_rel(r_ECI_SV1', v_ECI_SV1', r_ECI_SV3', v_ECI_SV3');
        SV3_RTN_pos = rho3';
        SV3_ECI_pos = SV3_state(1:3);

        % [rho3, ~] = ECI2RTN_rel(SV1_state(1:3)', SV1_state(4:6)', SV3_state(1:3)', SV3_state(4:6)');
        % SV3_RTN_pos = rho3';
        % SV3_ECI_pos = SV3_state(1:3);

        % TRY USING OTHER PROPAGATION FOR GROUND TRUTH
        
        RTN_noise = sqrtm(RTN_sigma)*randn(3,1); % use sigma value defined above
        ECI_noise = sqrtm(ECI_sigma)*randn(3,1);

        SV3_RTN_measurement = SV3_RTN_pos + RTN_noise;
        SV3_ECI_measurement = SV3_ECI_pos + ECI_noise;

        y_actual_EKF_SV3 = [SV3_RTN_measurement;SV3_ECI_measurement];

        % EKF Update
        H = measurement_sensitivity_matrix_SV3(SV1_state(1:3)',SV1_state(4:6)',SV1_OE_state);
        K = P_pred_EKF_SV3*(H')*inv(H*P_pred_EKF_SV3*H' + R);

        x_update_EKF_SV3_unscaled = x_pred_EKF_SV3_unscaled-(K*(y_actual_EKF_SV3-y_pred_EKF_SV3));
        x_update_EKF_SV3 = x_update_EKF_SV3_unscaled*a_chief;
        %x_update_EKF_SV3 = STM_curr*(x_update_EKF_SV3);
        %x_update_EKF_SV3 = ROEs_SV3_STM(i,:)';
        P_update_EKF_SV3 = (eye(6)-K*H)*P_pred_EKF_SV3*(eye(6)-K*H)' + K*R*K';

        % Extract updated chief QNS elements 
        a_qns_SV1     = SV1_OE_state(1);
        e_x_SV1       = SV1_OE_state(2);
        e_y_SV1       = SV1_OE_state(3);
        i_qns_SV1     = SV1_OE_state(4);
        RAAN_qns_SV1  = SV1_OE_state(5);
        u_SV1         = SV1_OE_state(6);
        
        % Extract deputy relative elements (from ROE state)
        d_a_SV3      = x_update_EKF_SV3(1);
        d_lambda_SV3 = x_update_EKF_SV3(2);
        d_e_x_SV3    = x_update_EKF_SV3(3);
        d_e_y_SV3    = x_update_EKF_SV3(4);
        d_i_x_SV3    = x_update_EKF_SV3(5);
        d_i_y_SV3    = x_update_EKF_SV3(6);
        
        % Reconstruct SV3's absolute QNS state 
        [a_qns_SV3, e_x_SV3, e_y_SV3, i_qns_SV3, RAAN_qns_SV3, u_SV3] = ...
            ROE2quasi_nonsing(a_qns_SV1, e_x_SV1, e_y_SV1, i_qns_SV1, RAAN_qns_SV1, u_SV1, ...
                              d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3);
        
        %SV3_OE_state = [a_qns_SV3, e_x_SV3, e_y_SV3, i_qns_SV3, RAAN_qns_SV3, u_SV3];

        % Save states to buffer
        state_SV3_all(i, :) = SV3_state';
        state_SV2_all(i, :) = SV2_state';
        state_SV1_all(i, :) = SV1_state';

        % EKF buffer
        OE_EKF_SV1_all(i, :) = SV1_OE_state';
        OE_EKF_SV3_all(i, :) = SV3_OE_state';
        x_pred_EKF_SV3_all(i, :) = x_pred_EKF_SV3;
        P_pred_EKF_SV3_all(i,:,:) = P_pred_EKF_SV3;
        y_actual_EKF_SV3_all(i,:) = y_actual_EKF_SV3;
        y_pred_EKF_SV3_all(i, :) = y_pred_EKF_SV3;
        x_update_EKF_SV3_all(i, :) = x_update_EKF_SV3;
        P_update_EKF_SV3_all(i,:,:) = P_update_EKF_SV3;
        

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
            %SV3_a_vals(i, :) = Lyapunov_feedback_control(SV3_roe, SV3_roe_nom, SV1_oe, N, k, num_orbits_station_keep(2)); % TODO: Fix this hardcoding
            % SV3_a_vals(i, :) = [0, 0, 0];
        end

        %%% APPLY ACCELERATION AS DELTA V 
        % Convert a_RTN to dv_RTN by multiplying by dt
        SV2_dv_vals(i,:) = dv_RTN2ECI(SV1_state(1:3), SV1_state(4:6), dt*SV2_a_vals(i, :)'/1e3); % m/s^2 --> m/s --> km/s
        SV3_dv_vals(i,:) = dv_RTN2ECI(SV1_state(1:3), SV1_state(4:6), dt*SV3_a_vals(i, :)'/1e3); % m/s^2 --> m/s --> km/s

        SV2_state(4:6) = SV2_state(4:6) + SV2_dv_vals(i,:)';
        %SV3_state(4:6) = SV3_state(4:6) + SV3_dv_vals(i,:)';

        delta_ROE_SV3 = delta_ROE_from_delta_v(dt*SV3_a_vals(i, :)',SV1_OE_state);
        %SV3_ROE_state = SV3_ROE_state + delta_ROE_SV3;
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
    % plot_3D_rel_orbit(SV2_rel_pos, SV3_rel_pos, true, 'All maneuvers 3D path', 'figures/PS7/ROE_3d_all_maneuvers.png');
    % % create_roe_3d_animation(SV2_rel_pos, SV3_rel_pos, 'SV2 and SV3 3D Relative Orbits (Animated)');
    % 
    RT = true; RN = true; NT = true;
    %plot_RT_RN_projections_separate(SV3_rel_pos, y_EKF_SV3_all(:,1:3), RT, RN, NT, '', 'figures/PS8/RTN_noise_comaprison.png');
    %plot_3D_orbits(r_SV3, y_actual_EKF_SV3_all(:,4:6), '', '', 'figures/PS8/ECI_noise_comparison.png')
    
    %[d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2] = ECI2ROE_array_mean(r_SV1, v_SV1, r_SV2, v_SV2, true); 
    [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3] = ECI2ROE_array(r_SV1, v_SV1, r_SV3, v_SV3); 

    ROE_SV3_true = [d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3];

    d_a_SV3_STM = ROEs_SV3_STM(:,1);
    d_lambda_SV3_STM =  ROEs_SV3_STM(:,2);
    d_e_x_SV3_STM = ROEs_SV3_STM(:,3);
    d_e_y_SV3_STM = ROEs_SV3_STM(:,4);
    d_i_x_SV3_STM = ROEs_SV3_STM(:,5);
    d_i_y_SV3_STM = ROEs_SV3_STM(:,6);

    d_a_SV3_EKF      = x_update_EKF_SV3_all(:, 1);
    d_lambda_SV3_EKF  = x_update_EKF_SV3_all(:, 2);
    d_e_x_SV3_EKF    = x_update_EKF_SV3_all(:, 3);
    d_e_y_SV3_EKF    = x_update_EKF_SV3_all(:, 4);
    d_i_x_SV3_EKF    = x_update_EKF_SV3_all(:, 5);
    d_i_y_SV3_EKF    = x_update_EKF_SV3_all(:, 6);

    plot_ROE_comparison(full_times, t_orbit, ...
    d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3, ...
    d_a_SV3_EKF, d_lambda_SV3_EKF, d_e_x_SV3_EKF, d_e_y_SV3_EKF, d_i_x_SV3_EKF, d_i_y_SV3_EKF,  ...
    'Ground Truth', 'EKF',  '','figures/PS8/ROE_planes_SV3_comparison.png', '', 'figures/PS8/ROE_over_time_SV3_comparison.png');
    plot_ROE_planes(full_times, t_orbit, d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3, 'figures/PS8/ROE_planes_ground_truth_SV3.png', 'figures/PS8/ROE_over_time_ground_truth_SV3.png');

    %plot_ROE_comparison(full_times, t_orbit, ...
    %d_a_SV3_STM, d_lambda_SV3_STM, d_e_x_SV3_STM, d_e_y_SV3_STM, d_i_x_SV3_STM, d_i_y_SV3_STM, ...
    %d_a_SV3_EKF, d_lambda_SV3_EKF, d_e_x_SV3_EKF, d_e_y_SV3_EKF, d_i_x_SV3_EKF, d_i_y_SV3_EKF,  ...
    %'Ground Truth', 'EKF Mean Propagation',  '', 'figures/PS8/ROE_planes_SV3_comparison.png', '', 'figures/PS8/ROE_over_time_SV3_comparison.png');
    plot_ROE_planes(full_times, t_orbit, d_a_SV3_STM, d_lambda_SV3_STM, d_e_x_SV3_STM, d_e_y_SV3_STM, d_i_x_SV3_STM, d_i_y_SV3_STM, 'figures/PS8/ROE_planes_ground_truth_STM_SV3.png', 'figures/PS8/ROE_over_time_ground_truth_STM_SV3.png');
    % %plot_ROE_planes(full_times, t_orbit, d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2, 'figures/PS6/ROE_planes_SV2.png', 'figures/PS6/ROE_over_time_SV2.png');
    % 
    % plot_ROE_planes_with_modes(full_times, t_orbit, d_a_SV2, d_lambda_SV2, d_e_x_SV2, d_e_y_SV2, d_i_x_SV2, d_i_y_SV2,SV2_modes,...
    %    num_orbits_modes, num_orbits_station_keep,'figures/PS7/ROE_planes_modes_SV2.png', 'figures/PS7/ROE_over_time_modes_SV2.png', 'figures/PS7/ROE_error_over_time_modes_SV2.png');
    % plot_ROE_planes_with_modes(full_times, t_orbit, d_a_SV3, d_lambda_SV3, d_e_x_SV3, d_e_y_SV3, d_i_x_SV3, d_i_y_SV3,SV3_modes,...
    %     num_orbits_modes, num_orbits_station_keep,'figures/PS7/ROE_planes_modes_SV3.png', 'figures/PS7/ROE_over_time_modes_SV3.png', 'figures/PS7/ROE_error_over_time_modes_SV3.png');
    % % plot_delta_v_timeline(full_times, SV2_a_vals, t_orbit, SV2_modes, num_orbits_modes, num_orbits_station_keep, ...
    %    true, 'figures/PS7/delta_v_timeline_modes_SV2.png');
    % plot_delta_v_timeline(full_times, SV3_a_vals, t_orbit, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
    %    true, 'figures/PS7/delta_v_timeline_modes_SV3.png');
    % plot_delta_v_cumulative_timeline(full_times, SV3_a_vals, t_orbit, SV3_modes, num_orbits_modes, num_orbits_station_keep, ...
    %    true, 'figures/PS7/delta_v_cumulative_timeline_modes_SV3.png')
    
    % figure_names = {
    % 'OE_EKF_SV1_all', 'OE_EKF_SV3_all', 'x_pred_EKF_SV3_all', ...
    % 'y_actual_EKF_SV3_all', 'y_pred_EKF_SV3_all', 'x_update_EKF_SV3_all'
    % };
    % 
    % data_vars = {
    %     OE_EKF_SV1_all, OE_EKF_SV3_all, x_pred_EKF_SV3_all, ...
    %     y_actual_EKF_SV3_all, y_pred_EKF_SV3_all, x_update_EKF_SV3_all
    % };
    % 
    % for idx = 1:length(data_vars)
    %     data = data_vars{idx};
    %     name = figure_names{idx};
    % 
    %     figure('Name', name, 'NumberTitle', 'off');
    %     for k = 1:6
    %         subplot(3,2,k);
    %         plot(full_times, data(:,k), 'LineWidth', 1.2);
    %         xlabel('Time Step');
    %         ylabel(sprintf('State %d', k));
    %         title(sprintf('%s - Element %d', name, k));
    %         grid on;
    %     end
    %     sgtitle(strrep(name, '_', '\_'));  % Add overall title
    % end

end