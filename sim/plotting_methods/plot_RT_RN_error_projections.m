function plot_RT_RN_error_projections(SV2_rel_pos1, SV2_rel_pos2, SV2_rel_pos_true, SV3_rel_pos1, SV3_rel_pos2, SV3_rel_pos_true, title_string,fig_path)
    % Extract coordinates
    x1_true = SV2_rel_pos_true(:,1); y1_true = SV2_rel_pos_true(:,2); z1_true = SV2_rel_pos_true(:,3);
    x2_true = SV3_rel_pos_true(:,1); y2_true = SV3_rel_pos_true(:,2); z2_true = SV3_rel_pos_true(:,3);
    x1_1 = SV2_rel_pos1(:,1)-x1_true; y1_1 = SV2_rel_pos1(:,2)-y1_true; z1_1 = SV2_rel_pos1(:,3)-z1_true;
    x2_1 = SV3_rel_pos1(:,1)-x2_true; y2_1 = SV3_rel_pos1(:,2)-y2_true; z2_1 = SV3_rel_pos1(:,3)-z2_true;
    x1_2 = SV2_rel_pos2(:,1)-x1_true; y1_2 = SV2_rel_pos2(:,2)-y1_true; z1_2 = SV2_rel_pos2(:,3)-z1_true;
    x2_2 = SV3_rel_pos2(:,1)-x2_true; y2_2 = SV3_rel_pos2(:,2)-y2_true; z2_2 = SV3_rel_pos2(:,3)-z2_true;

    figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
    
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % === X-Y Plane ===
    ax1 = nexttile;
    plot(y1_1, x1_1, 'c', 'LineWidth', 1.5); hold on;
    plot(y1_2, x1_2, 'g', 'LineWidth', 1.5); 
    plot(y2_1, x2_1, 'm--', 'LineWidth', 1.5);
    plot(y2_2, x2_2, 'y--', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    axis equal;
    ylabel('Radial Error [km]'); xlabel('Tangential Error [km]');
    title('R-T Plane Projection');
    
    % === X-Z Plane ===
    ax2 = nexttile;
    plot(z1_1, x1_1, 'c', 'LineWidth', 1.5); hold on;
    plot(z1_2, x1_2, 'g', 'LineWidth', 1.5); 
    plot(z2_1, x2_1, 'm--', 'LineWidth', 1.5);
    plot(z2_2, x2_2, 'y--', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    axis equal;
    ylabel('Radial Error [km]'); xlabel('Normal Error [km]');
    title('R-N Plane Projection');
    
    % === Y-Z Plane ===
    ax3 = nexttile;
    plot(y1_1, z1_1, 'c', 'LineWidth', 1.5); hold on;
    plot(y1_2, z1_2, 'g', 'LineWidth', 1.5); 
    plot(y2_1, z2_1, 'm--', 'LineWidth', 1.5);
    plot(y2_2, z2_2, 'y--', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    axis equal;
    ylabel('Normal Error [km]'); xlabel('Tangential Error [km]');
    title('N-T Plane Projection');
    
    % === Global Legend ===
    lgd = legend({'SV2 / Watcher with YA Dif. Eq.', 'SV2 / Watcher with Geometric Mapping', ...
        'SV3 / Watcher with YA Dif. Eq.', 'SV3 / Watcher with Geometric Mapping', ...
        'SV1 / Chief'}, 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 3);
    lgd.Layout.Tile = 'south'; % Place legend below all subplots
    
    % Save figure
    saveas(gcf, fig_path);
end
