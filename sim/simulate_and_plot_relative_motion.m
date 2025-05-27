function simulate_and_plot_relative_motion(t, t_orbit, SV2_roe_init, state_abs_SV2_init, SV3_roe_init, state_abs_SV3_init, SV1_oe_init, SV1_state_init,...
    fig_path, title_str, fig_path2, title_str2, fig_path3, title_str3)
    % Inputs:
    %   t             = time array [s]
    %   SV2_roe_init  = 6x1 initial QNS ROEs of Watcher (SV2)
    %   SV3_roe_init  = 6x1 initial QNS ROEs of Docker  (SV3)
    %   SV1_oe_init   = 6x1 initial classical orbital elements of Chief (a, e, i, Omega, w, M) in degrees
    %   fig_path      = full path to save figure (e.g., 'figs/formation_projection.png')
    %   title_str     = string to display in the plot title

    J2_flag = true;

    % Propagate chief
    [~, chief_state] = rk4_eom_ECI(t, SV1_state_init, J2_flag); % with J2
    r_SV1 = chief_state(:,1:3);
    v_SV1 = chief_state(:,4:6);
    
    % Propagate both deputies in ECI
    [~, state_SV2] = rk4_eom_ECI(t, state_abs_SV2_init, J2_flag); % with J2
    [~, state_SV3] = rk4_eom_ECI(t, state_abs_SV3_init, J2_flag); % with J2
    r_SV2 = state_SV2(:,1:3);
    v_SV2 = state_SV2(:,4:6);
    r_SV3 = state_SV3(:,1:3);
    v_SV3 = state_SV3(:,4:6);

    % RTN positions from FODE
    SV2_rel_pos = zeros(length(t), 3);
    SV3_rel_pos = zeros(length(t), 3);

    % ROE from FODE
    [d_a_osc_SV2, d_lambda_osc_SV2, d_e_x_osc_SV2, d_e_y_osc_SV2, d_i_x_osc_SV2, d_i_y_osc_SV2, ...
     d_a_mean_SV2, d_lambda_mean_SV2, d_e_x_mean_SV2, d_e_y_mean_SV2, d_i_x_mean_SV2, d_i_y_mean_SV2] = ...
        compute_OE_ROE_mean_osc(r_SV1, v_SV1, r_SV2, v_SV2, J2_flag);
    [d_a_osc_SV3, d_lambda_osc_SV3, d_e_x_osc_SV3, d_e_y_osc_SV3, d_i_x_osc_SV3, d_i_y_osc_SV3, ...
     d_a_mean_SV3, d_lambda_mean_SV3, d_e_x_mean_SV3, d_e_y_mean_SV3, d_i_x_mean_SV3, d_i_y_mean_SV3] = ...
        compute_OE_ROE_mean_osc(r_SV1, v_SV1, r_SV3, v_SV3, J2_flag);
    ROE_SV2_FODE = [d_a_mean_SV2, d_lambda_mean_SV2, d_e_x_mean_SV2, d_e_y_mean_SV2, d_i_x_mean_SV2, d_i_y_mean_SV2];
    ROE_SV3_FODE = [d_a_mean_SV3, d_lambda_mean_SV3, d_e_x_mean_SV3, d_e_y_mean_SV3, d_i_x_mean_SV3, d_i_y_mean_SV3];

    % Propagate ROE from STM using same IC

    [~, roe_SV2] = roe_stm_j2(t, ROE_SV2_FODE(1,:), SV1_oe_init);
    [~, roe_SV3] = roe_stm_j2(t, ROE_SV3_FODE(1,:), SV1_oe_init);


    for i = 1:length(t)     
        % Unpack ECI states for each sc at each timestep
        x_r_SV1 = r_SV1(i,:)';
        x_v_SV1 = v_SV1(i,:)';
        x_r_SV2 = r_SV2(i,:)';
        x_v_SV2 = v_SV2(i,:)';
        x_r_SV3 = r_SV3(i,:)';
        x_v_SV3 = v_SV3(i,:)';

        SV2_rel_pos(i, :) = ECI2RTN_rel(x_r_SV1', x_v_SV1', x_r_SV2', x_v_SV2');
        SV3_rel_pos(i, :) = ECI2RTN_rel(x_r_SV1', x_v_SV1', x_r_SV3', x_v_SV3');     

    end

    % Plot 3-view projections
    RT = true;
    RN = true;
    NT = true; 
    plot_RT_RN_projections_separate(SV2_rel_pos, SV3_rel_pos, RT,RN,NT,title_str, fig_path);
    %plot_RT_RN_projections(SV2_rel_pos, SV3_rel_pos,title_str, fig_path);

    plot_ROE_comparison(t, t_orbit, roe_SV2, ROE_SV2_FODE, 'STM', 'Dif. FODE', title_str2, fig_path2, title_str3, fig_path3)
end