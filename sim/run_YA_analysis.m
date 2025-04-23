function [SV2_rel_pos_ecc,SV2_rel_vel_ecc,SV3_rel_pos_ecc,SV3_rel_vel_ecc,...
    SV2_YA_pos,SV2_YA_vel,SV3_YA_pos,SV3_YA_vel,SV2_YA_mapping_pos,SV2_YA_mapping_vel,SV3_YA_mapping_pos,SV3_YA_mapping_vel,...
    rho_pos_ratio_SV2_avg,rho_pos_ratio_SV3_avg,K_YA_SV2,K_YA_SV3,ROE_SV2_unscaled,ROE_SV3_unscaled] = ...
    run_YA_analysis(tstart,tint,tend,r_RTN_no_j2,a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init,e_SV1_init,w_SV1_init,nu_SV1_init,M_SV1_init,r_SV1_ECI_init, v_SV1_ECI_init, ...
    d_a_SV2_init,d_lambda_SV2_init,d_e_x_SV2_init,d_e_y_SV2_init,d_i_x_SV2_init,d_i_y_SV2_init,...
    d_a_SV3_init,d_lambda_SV3_init,d_e_x_SV3_init,d_e_y_SV3_init,d_i_x_SV3_init,d_i_y_SV3_init)
    
    [r_SV2_init, v_SV2_init] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
    d_a_SV2_init,d_lambda_SV2_init,d_e_x_SV2_init,d_e_y_SV2_init,d_i_x_SV2_init,d_i_y_SV2_init);
    state_abs_SV2_init = [r_SV2_init;v_SV2_init];
    
    [r_SV3_init, v_SV3_init] = ROE2ECI(a_SV1_init,ex_SV1_init,ey_SV1_init,i_SV1_init,RAAN_SV1_init,u_SV1_init, ...
    d_a_SV3_init,d_lambda_SV3_init,d_e_x_SV3_init,d_e_y_SV3_init,d_i_x_SV3_init,d_i_y_SV3_init);
    state_abs_SV3_init = [r_SV3_init;v_SV3_init];
    
    % Converting initial swarm ECI coordinates to initial relative co-ordinates in
    % chief's RTN frame.
    [r_SV2_RTN_init, v_SV2_RTN_init] = ECI2RTN_rel(r_SV1_ECI_init', v_SV1_ECI_init', r_SV2_init', v_SV2_init');
    r_SV2_RTN_init = r_SV2_RTN_init';
    v_SV2_RTN_init = v_SV2_RTN_init';
    
    [r_SV3_RTN_init, v_SV3_RTN_init] = ECI2RTN_rel(r_SV1_ECI_init', v_SV1_ECI_init', r_SV3_init', v_SV3_init');
    r_SV3_RTN_init = r_SV3_RTN_init';
    v_SV3_RTN_init = v_SV3_RTN_init';
    
    % This needs to also have the chief's state to propagate
    state_rel_SV2_init = [r_SV1_ECI_init; v_SV1_ECI_init; r_SV2_RTN_init; v_SV2_RTN_init];
    state_rel_SV3_init = [r_SV1_ECI_init; v_SV1_ECI_init; r_SV3_RTN_init; v_SV3_RTN_init];
    
    %%%%%% RUN SIM OF RELATIVE MOTION %%%%%%%%
    [t_2, state_2] = rk4_eom_rel_RTN(tstart:tint:tend, state_rel_SV2_init);
    [t_3, state_3] = rk4_eom_rel_RTN(tstart:tint:tend, state_rel_SV3_init);
    SV2_rel_pos_ecc = [state_2(:, 7), state_2(:, 8), state_2(:, 9)];
    SV2_rel_vel_ecc = [state_2(:, 10), state_2(:, 11), state_2(:, 12)];
    SV3_rel_pos_ecc = [state_3(:, 7), state_3(:, 8), state_3(:, 9)];
    SV3_rel_vel_ecc = [state_3(:, 10), state_3(:, 11), state_3(:, 12)];
    
    % We see that the ratio of relative distance to orbit size is smaller than 1e-3
    rho_pos_ratio_SV2_avg = mean(vecnorm(SV2_rel_pos_ecc, 2, 2)./r_RTN_no_j2(:, 1));
    rho_pos_ratio_SV3_avg = mean(vecnorm(SV3_rel_pos_ecc, 2, 2)./r_RTN_no_j2(:, 1));
    
    
    %%%%% CALCULATE THE YA INTEGRATION CONSTANTS %%%%% 
    K_YA_SV2 = ya_integration_constants([r_SV2_RTN_init; v_SV2_RTN_init], a_SV1_init, e_SV1_init, nu_SV1_init);
    K_YA_SV3 = ya_integration_constants([r_SV3_RTN_init; v_SV3_RTN_init], a_SV1_init, e_SV1_init, nu_SV1_init);
    
    
    %%%% PROPAGATE RELATIVE POS, VEL USING YA SOLUTION %%%%%%
    % evaluate_ya(t, a, e, K, M_init)
    SV2_YA_state = evaluate_ya(t_2, a_SV1_init, e_SV1_init, M_SV1_init, K_YA_SV2);
    SV3_YA_state = evaluate_ya(t_2, a_SV1_init, e_SV1_init, M_SV1_init, K_YA_SV3);  
    
    SV2_YA_pos = SV2_YA_state(:, 1:3);
    SV3_YA_pos = SV3_YA_state(:, 1:3);
    SV2_YA_vel = SV2_YA_state(:, 4:6);
    SV3_YA_vel = SV3_YA_state(:, 4:6);
    
    %%%% PROPAGATE RELATIVE POS, VEL USING YA GEOMETRIC LINEAR MAPPING %%%%%%
    % Assume ROE are passed in as meters, need to divide by chief's
    % semi-major axis in meters to get into expression for geometric mapping

    a_chief = a_SV1_init*1e3; % km --> m

    d_a_SV2 = d_a_SV2_init/a_chief; d_lambda_SV2 = d_lambda_SV2_init/a_chief; d_e_x_SV2 = d_e_x_SV2_init/a_chief;
    d_e_y_SV2 = d_e_y_SV2_init/a_chief; d_i_x_SV2 = d_i_x_SV2_init/a_chief; d_i_y_SV2 = d_i_y_SV2_init/a_chief;

    ROE_SV2_unscaled = [d_a_SV2;d_lambda_SV2;d_e_x_SV2;d_e_y_SV2;d_i_x_SV2;d_i_y_SV2];

    d_a_SV3 = d_a_SV3_init/a_chief; d_lambda_SV3 = d_lambda_SV3_init/a_chief; d_e_x_SV3 = d_e_x_SV3_init/a_chief;
    d_e_y_SV3 = d_e_y_SV3_init/a_chief; d_i_x_SV3 = d_i_x_SV3_init/a_chief; d_i_y_SV3 = d_i_y_SV3_init/a_chief;

    ROE_SV3_unscaled = [d_a_SV3;d_lambda_SV3;d_e_x_SV3;d_e_y_SV3;d_i_x_SV3;d_i_y_SV3];

    SV2_YA_state_mapping = evaluate_YA_geometric_mapping(t_2,a_SV1_init,ex_SV1_init,...
        ey_SV1_init,i_SV1_init,w_SV1_init,M_SV1_init,d_a_SV2,d_lambda_SV2,...
        d_e_x_SV2,d_e_y_SV2,d_i_x_SV2,d_i_y_SV2);
    SV3_YA_state_mapping = evaluate_YA_geometric_mapping(t_2,a_SV1_init,ex_SV1_init,...
        ey_SV1_init,i_SV1_init,w_SV1_init,M_SV1_init,d_a_SV3,d_lambda_SV3,...
        d_e_x_SV3,d_e_y_SV3,d_i_x_SV3,d_i_y_SV3);
    
    SV2_YA_mapping_pos = SV2_YA_state_mapping(:, 1:3);
    SV3_YA_mapping_pos = SV3_YA_state_mapping(:, 1:3);
    SV2_YA_mapping_vel = SV2_YA_state_mapping(:, 4:6);
    SV3_YA_mapping_vel = SV3_YA_state_mapping(:, 4:6);
end