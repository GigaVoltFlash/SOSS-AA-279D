function plot_RTN_error_vs_time(SV2_rel_pos1, SV2_rel_pos2, SV2_rel_pos_true, ...
                                      SV3_rel_pos1, SV3_rel_pos2, SV3_rel_pos_true, ...
                                      time_vec,  fig_path_SV2, ...
                                      fig_path_SV3)
% Plot R,T,N error over time for both methods, separately for SV2 and SV3

    % --- Compute errors for SV2 ---
    err_R_1_SV2 = SV2_rel_pos1(:,1) - SV2_rel_pos_true(:,1);
    err_T_1_SV2 = SV2_rel_pos1(:,2) - SV2_rel_pos_true(:,2);
    err_N_1_SV2 = SV2_rel_pos1(:,3) - SV2_rel_pos_true(:,3);

    err_R_2_SV2 = SV2_rel_pos2(:,1) - SV2_rel_pos_true(:,1);
    err_T_2_SV2 = SV2_rel_pos2(:,2) - SV2_rel_pos_true(:,2);
    err_N_2_SV2 = SV2_rel_pos2(:,3) - SV2_rel_pos_true(:,3);

    % --- Compute errors for SV3 ---
    err_R_1_SV3 = SV3_rel_pos1(:,1) - SV3_rel_pos_true(:,1);
    err_T_1_SV3 = SV3_rel_pos1(:,2) - SV3_rel_pos_true(:,2);
    err_N_1_SV3 = SV3_rel_pos1(:,3) - SV3_rel_pos_true(:,3);

    err_R_2_SV3 = SV3_rel_pos2(:,1) - SV3_rel_pos_true(:,1);
    err_T_2_SV3 = SV3_rel_pos2(:,2) - SV3_rel_pos_true(:,2);
    err_N_2_SV3 = SV3_rel_pos2(:,3) - SV3_rel_pos_true(:,3);

    % ================================================================
    % Figure 1 — SV2
    % ================================================================
    figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    %sgtitle(title_string_SV2, 'FontWeight', 'bold');

    % R
    ax1 = nexttile;
    plot(time_vec, err_R_1_SV2, 'c', 'LineWidth', 1.5); hold on;
    plot(time_vec, err_R_2_SV2, 'g', 'LineWidth', 1.5);
    grid on;
    xlabel('Time'); ylabel('Radial Error [km]');
    title('Radial (R) Error vs Time');

    % T
    ax2 = nexttile;
    plot(time_vec, err_T_1_SV2, 'c', 'LineWidth', 1.5); hold on;
    plot(time_vec, err_T_2_SV2, 'g', 'LineWidth', 1.5);
    grid on;
    xlabel('Time'); ylabel('Tangential Error [km]');
    title('Tangential (T) Error vs Time');

    % N
    ax3 = nexttile;
    plot(time_vec, err_N_1_SV2, 'c', 'LineWidth', 1.5); hold on;
    plot(time_vec, err_N_2_SV2, 'g', 'LineWidth', 1.5);
    grid on;
    xlabel('Time'); ylabel('Normal Error [km]');
    title('Normal (N) Error vs Time');

    % Legend
    lgd = legend({'YA Dif. Eq.', 'Geometric Mapping'}, ...
                  'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2);
    lgd.Layout.Tile = 'south';

    % Save
    saveas(gcf, fig_path_SV2);

    % ================================================================
    % Figure 2 — SV3
    % ================================================================
    figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    %sgtitle(title_string_SV3, 'FontWeight', 'bold');

    % R
    ax1 = nexttile;
    plot(time_vec, err_R_1_SV3, 'm', 'LineWidth', 1.5); hold on;
    plot(time_vec, err_R_2_SV3, 'y', 'LineWidth', 1.5);
    grid on;
    xlabel('Time'); ylabel('Radial Error [km]');
    title('Radial (R) Error vs Time');

    % T
    ax2 = nexttile;
    plot(time_vec, err_T_1_SV3, 'm', 'LineWidth', 1.5); hold on;
    plot(time_vec, err_T_2_SV3, 'y', 'LineWidth', 1.5);
    grid on;
    xlabel('Time'); ylabel('Tangential Error [km]');
    title('Tangential (T) Error vs Time');

    % N
    ax3 = nexttile;
    plot(time_vec, err_N_1_SV3, 'm', 'LineWidth', 1.5); hold on;
    plot(time_vec, err_N_2_SV3, 'y', 'LineWidth', 1.5);
    grid on;
    xlabel('Time'); ylabel('Normal Error [km]');
    title('Normal (N) Error vs Time');

    % Legend
    lgd = legend({'YA Dif. Eq.', 'Geometric Mapping'}, ...
                  'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2);
    lgd.Layout.Tile = 'south';

    % Save
    saveas(gcf, fig_path_SV3);
end
