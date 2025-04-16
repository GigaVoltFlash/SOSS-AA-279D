function plot_rel_sat_pos_3D(SV2_rel_pos, SV3_rel_pos, fig_path)
    % Inputs:
    % SV2_rel_pos - nx3 array of [x, y, z] for SV2/Watcher
    % SV3_rel_pos - nx3 array of [x, y, z] for SV3/Docker

    % Extract coordinates
    x1 = SV2_rel_pos(:,1); y1 = SV2_rel_pos(:,2); z1 = SV2_rel_pos(:,3);
    x2 = SV3_rel_pos(:,1); y2 = SV3_rel_pos(:,2); z2 = SV3_rel_pos(:,3);

    % Create figure
    figure('Color', 'w');
    
    % SV2 / Watcher subplot
    subplot(1,2,1);
    plot3(x1, y1, z1, 'b', 'LineWidth', 1.5); hold on;
    % scatter3(0, 0, 0, 50, 'ko', 'filled'); % Mark origin
    text(x1(end), y1(end), z1(end), ' SV2/Watcher', 'Color', 'b', 'FontSize', 10);
    % text(0, 0, 0, ' SV1/Target/Chief', 'FontSize', 9, 'Color', 'k');
    grid on; axis equal; view(3);
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title('SV2/Watcher Relative Trajectory');

    % SV3 / Docker subplot
    subplot(1,2,2);
    plot3(x2, y2, z2, 'r', 'LineWidth', 1.5); hold on;
    % scatter3(0, 0, 0, 50, 'ko', 'filled'); % Mark origin
    text(x2(end), y2(end), z2(end), ' SV3/Docker', 'Color', 'r', 'FontSize', 10);
    % text(0, 0, 0, ' SV1/Target/Chief', 'FontSize', 9, 'Color', 'k');
    grid on; axis equal; view(3);
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title('SV3/Docker Relative Trajectory');

    % Super title
    sgtitle('Relative Positions of SV2 and SV3 in RTN Frame');

    % Save figure
    saveas(gcf, fig_path);
end
