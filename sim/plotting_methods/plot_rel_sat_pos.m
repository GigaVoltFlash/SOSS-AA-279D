function plot_rel_sat_pos(SV2_rel_pos, SV3_rel_pos)
    % Inputs:
    % SV2_rel_pos - nx3 array of [x, y, z] for SV2/Watcher
    % SV3_rel_pos - nx3 array of [x, y, z] for SV3/Docker

    % Extract coordinates
    x1 = SV2_rel_pos(:,1);
    y1 = SV2_rel_pos(:,2);
    z1 = SV2_rel_pos(:,3);

    x2 = SV3_rel_pos(:,1);
    y2 = SV3_rel_pos(:,2);
    z2 = SV3_rel_pos(:,3);

    % Plotting
    figure('Color', 'w');
    hold on; grid on; axis equal;

    % Plot trajectories
    plot3(x1, y1, z1, 'b', 'LineWidth', 1.5);
    plot3(x2, y2, z2, 'r', 'LineWidth', 1.5);

    % Final point labels
    text(x1(end), y1(end), z1(end), ' SV2/Watcher', 'Color', 'b', 'FontSize', 10);
    text(x2(end), y2(end), z2(end), ' SV3/Docker', 'Color', 'r', 'FontSize', 10);

    % Mark and label the origin
    scatter3(0, 0, 0, 70, 'ko', 'filled');
    text(0, 0, 0, ' SV1/Target/Chief', 'Color', 'k', 'FontSize', 10, 'HorizontalAlignment', 'left');

    % Labels and title
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    title('Relative Satellite Positions in 3D (RTN Frame)');
    legend({'SV2/Watcher', 'SV3/Docker'}, 'Location', 'best');

    % Save the figure
    if ~exist('figures', 'dir')
        mkdir('figures');
    end
    saveas(gcf, '../figures/satellite_trajectories.png');
end
