function simulate_and_plot_relative_motion(t, SV2_roe_init, SV3_roe_init, SV1_oe_init, fig_path, title_str)
    % Inputs:
    %   t             = time array [s]
    %   SV2_roe_init  = 6x1 initial QNS ROEs of Watcher (SV2)
    %   SV3_roe_init  = 6x1 initial QNS ROEs of Docker  (SV3)
    %   SV1_oe_init   = 6x1 classical orbital elements of Chief (a, e, i, Omega, w, M) in degrees
    %   fig_path      = full path to save figure (e.g., 'figs/formation_projection.png')
    %   title_str     = string to display in the plot title

    global mu_earth;

    % Propagate both deputies
    [~, roe_SV2] = roe_stm_j2(t, SV2_roe_init, SV1_oe_init);
    [~, roe_SV3] = roe_stm_j2(t, SV3_roe_init, SV1_oe_init);

    % Chief parameters
    a = SV1_oe_init(1);     % semi-major axis [m]
    n = sqrt(mu_earth / a^3);

    % Convert ROEs to RTN coordinates (in kilometers)
    SV2_rel_pos = zeros(length(t), 3);
    SV3_rel_pos = zeros(length(t), 3);
    
    for i = 1:length(t)
        % Unpack ROEs for each deputy
        x_roe_SV2 = roe_SV2(i, :)';
        x_roe_SV3 = roe_SV3(i, :)';

        % Convert QNS ROEs to relative RTN position
        SV2_rel_pos(i, :) = ROE2RTN(x_roe_SV2, SV1_oe_init);
        SV3_rel_pos(i, :) = ROE2RTN(x_roe_SV3, SV1_oe_init);
    end

    % Plot 3-view projections
    plot_RT_RN_projections(SV2_rel_pos, SV3_rel_pos, title_str, fig_path);
end