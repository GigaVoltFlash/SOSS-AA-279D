function plot_RT_RN_projections_both_PS2_2(SV2_rel_pos1, SV3_rel_pos1, title_string, fig_path)
    % Extract coordinates
    x1_1 = SV2_rel_pos1(:,1); y1_1 = SV2_rel_pos1(:,2); z1_1 = SV2_rel_pos1(:,3);
    x2_1 = SV3_rel_pos1(:,1); y2_1 = SV3_rel_pos1(:,2); z2_1 = SV3_rel_pos1(:,3);

    % === Setup tiledlayout ===
    figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle(title_string, 'FontWeight', 'bold');

    % === X-Y Plane ===
    ax1 = nexttile;
    p1 = plot(y1_1, x1_1, 'g', 'LineWidth', 3); hold on;
    %p2 = plot(y1_2, x1_2, 'b--', 'LineWidth', 1.5); 
    p3 = plot(y2_1, x2_1, 'm', 'LineWidth', 3);
    %p4 = plot(y2_2, x2_2, 'r--', 'LineWidth', 1.5);
    p5 = plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    ylabel('Radial [km]'); xlabel('Tangential [km]');
    title('R-T Plane Projection');

    % === X-Z Plane ===
    ax2 = nexttile;
    plot(z1_1, x1_1, 'g', 'LineWidth', 3); hold on;
    %plot(z1_2, x1_2, 'b--', 'LineWidth', 1.5); 
    plot(z2_1, x2_1, 'm', 'LineWidth', 3);
    %plot(z2_2, x2_2, 'r--', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    ylabel('Radial [km]'); xlabel('Normal [km]');
    title('R-N Plane Projection');

    % === Y-Z Plane ===
    ax3 = nexttile;
    plot(y1_1, z1_1, 'g', 'LineWidth', 3); hold on;
    %plot(y1_2, z1_2, 'b--', 'LineWidth', 1.5); 
    plot(y2_1, z2_1, 'm', 'LineWidth', 3);
    %plot(y2_2, z2_2, 'r--', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    ylabel('Normal [km]'); xlabel('Tangential [km]');
    title('N-T Plane Projection');

    % === Global Legend ===
    lgd = legend([p1,  p3,  p5], ...
        {'SV2 / Watcher with Maneuver',  ...
         'SV3 / Docker with with Maneuver', 'SV1 / Chief'}, ...
         'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 3);
    lgd.Layout.Tile = 'south';

    % Save figure
    saveas(gcf, fig_path);
end
