function plot_rel_sat_pos_3D_multi_method(SV2_rel_pos1, SV2_rel_pos2, SV3_rel_pos1, SV3_rel_pos2, fig_path)

    x1_1 = SV2_rel_pos1(:,1); y1_1 = SV2_rel_pos1(:,2); z1_1 = SV2_rel_pos1(:,3);
    x2_1 = SV3_rel_pos1(:,1); y2_1 = SV3_rel_pos1(:,2); z2_1 = SV3_rel_pos1(:,3);
    x1_2 = SV2_rel_pos2(:,1); y1_2 = SV2_rel_pos2(:,2); z1_2 = SV2_rel_pos2(:,3);
    x2_2 = SV3_rel_pos2(:,1); y2_2 = SV3_rel_pos2(:,2); z2_2 = SV3_rel_pos2(:,3);

    % Create figure
    figure('Color', 'w', 'Position', [100, 100, 1200, 600]);
    
    % SV2 / Watcher subplot
    subplot(1,2,1);
    plot3(x1_1, y1_1, z1_1, 'b', 'LineWidth', 1.5); hold on;
    plot3(x1_2, y1_2, z1_2, 'r--', 'LineWidth', 1.5); hold on;
    % scatter3(0, 0, 0, 50, 'ko', 'filled'); % Mark origin
    text(x1_1(end), y1_1(end), z1_1(end), ' SV2/Watcher with relative EOM', 'Color', 'b', 'FontSize', 10);
    text(x1_2(end-100), y1_2(end-100), z1_2(end-100), ' SV2/Watcher with FODE', 'Color', 'r', 'FontSize', 10);
    % text(0, 0, 0, ' SV1/Target/Chief', 'FontSize', 9, 'Color', 'k');
    grid on; axis equal; view(3);
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

    % SV3 / Docker subplot
    subplot(1,2,2);
    plot3(x2_1, y2_1, z2_1, 'b', 'LineWidth', 1.5); hold on;
    plot3(x2_2, y2_2, z2_2, 'r--', 'LineWidth', 1.5); hold on;
    % scatter3(0, 0, 0, 50, 'ko', 'filled'); % Mark origin
    text(x2_1(end), y2_1(end), z2_1(end), ' SV3/Docker with relative EOM', 'Color', 'b', 'FontSize', 10);
    text(x2_2(end-100), y2_2(end-100), z2_2(end-100), ' SV3/Docker with FODE', 'Color', 'r', 'FontSize', 10);
    % text(0, 0, 0, ' SV1/Target/Chief', 'FontSize', 9, 'Color', 'k');
    grid on; axis equal; view(3);
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

    % Super title
    sgtitle('Relative Positions of SV2 and SV3 in RTN Frame');

    % Save figure
    saveas(gcf, fig_path);
end
