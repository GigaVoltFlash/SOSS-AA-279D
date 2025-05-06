function simulate_and_plot_relative_motion(t, SV2_roe_init, state_abs_SV2_init, SV3_roe_init, state_abs_SV3_init, SV1_oe_init, SV1_state_init, fig_path, title_str)
    % Inputs:
    %   t             = time array [s]
    %   SV2_roe_init  = 6x1 initial QNS ROEs of Watcher (SV2)
    %   SV3_roe_init  = 6x1 initial QNS ROEs of Docker  (SV3)
    %   SV1_oe_init   = 6x1 initial classical orbital elements of Chief (a, e, i, Omega, w, M) in degrees
    %   fig_path      = full path to save figure (e.g., 'figs/formation_projection.png')
    %   title_str     = string to display in the plot title

    % Propagate chief
    [~, chief_state] = rk4_eom_ECI(t, SV1_state_init, true); % with J2
    r_SV1 = chief_state(:,1:3);
    v_SV1 = chief_state(:,4:6);
    
    % Propagate both deputies 
    % NEED TO FIX
    [~, roe_SV2] = roe_stm_j2(t, SV2_roe_init, SV1_oe_init);
    [~, roe_SV3] = roe_stm_j2(t, SV3_roe_init, SV1_oe_init);

    [~, state_SV2] = rk4_eom_ECI(t, state_abs_SV2_init, true); % with J2
    [~, state_SV3] = rk4_eom_ECI(t, state_abs_SV3_init, true); % with J2
    r_SV2 = state_SV2(:,1:3);
    v_SV2 = state_SV2(:,4:6);
    r_SV3 = state_SV3(:,1:3);
    v_SV3 = state_SV3(:,4:6);

    % Convert ROEs to RTN coordinates (in kilometers)
    SV2_rel_pos = zeros(length(t), 3);
    SV3_rel_pos = zeros(length(t), 3);
    
    for i = 1:length(t)
        % Unpack ROEs for each deputy
        x_roe_SV2 = roe_SV2(i, :)';
        x_roe_SV3 = roe_SV3(i, :)';

        x_r_SV1 = r_SV1(i,:)';
        x_v_SV1 = v_SV1(i,:)';
        x_r_SV2 = r_SV2(i,:)';
        x_v_SV2 = v_SV2(i,:)';
        x_r_SV3 = r_SV3(i,:)';
        x_v_SV3 = v_SV3(i,:)';

        [a_o, e_o, i_o, RAAN_o, w_o, M_o] = ECI2OE(x_r_SV1, x_v_SV1);
        oe_SV1 = [a_o, e_o, i_o, RAAN_o, w_o, M_o];

        % Convert QNS ROEs to relative RTN position
        %SV2_rel_pos(i, :) = ROE2RTN(x_roe_SV2, oe_SV1,x_r_SV1,x_v_SV1);
        %SV3_rel_pos(i, :) = ROE2RTN(x_roe_SV3, oe_SV1,x_r_SV1,x_v_SV1);
        % For now, just propagating regularly, not using J2 STM
        SV2_rel_pos(i, :) = ECI2RTN_rel(x_r_SV1', x_v_SV1', x_r_SV2', x_v_SV2');
        SV3_rel_pos(i, :) = ECI2RTN_rel(x_r_SV1', x_v_SV1', x_r_SV3', x_v_SV3');


    end

    % Plot 3-view projections
    RT = true;
    RN = true;
    NT = true; 
    plot_RT_RN_projections_separate(SV2_rel_pos, SV3_rel_pos, RT,RN,NT,title_str, fig_path);
    %plot_RT_RN_projections(SV2_rel_pos, SV3_rel_pos,title_str, fig_path);
end