function plot_RT_RN_projections(SV2_rel_pos, SV3_rel_pos, title_string,fig_path)
    % Extract coordinates
    x1 = SV2_rel_pos(:,1); y1 = SV2_rel_pos(:,2); z1 = SV2_rel_pos(:,3);
    x2 = SV3_rel_pos(:,1); y2 = SV3_rel_pos(:,2); z2 = SV3_rel_pos(:,3);

    figure('Color', 'w');
    
    % === X-Y Plane ===
    subplot(1, 3, 1);
    plot(y1, x1, 'b', 'LineWidth', 1.5); hold on;
    plot(y2, x2, 'r', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    axis equal;
    ylabel('Radial [km]'); xlabel('Tangential [km]');
    title('R-T Plane Projection');
    legend('SV2 / Watcher', 'SV3 / Docker', 'SV1 / Chief');
    % legend('SV2 / Watcher', 'SV1 / Chief');

    % === X-Z Plane ===
    subplot(1, 3, 2);
    plot(z1, x1, 'b', 'LineWidth', 1.5); hold on;
    plot(z2, x2, 'r', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    axis equal;
    ylabel('Radial [km]'); xlabel('Normal [km]');
    title('R-N Plane Projection');
    legend('SV2 / Watcher', 'SV3 / Docker', 'SV1 / Chief');

    % === Y-Z Plane ===
    subplot(1, 3, 3);
    plot(y1, z1, 'b', 'LineWidth', 1.5); hold on;
    plot(y2, z2, 'r', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    axis equal;
    ylabel('Normal [km]'); xlabel('Tangential [km]');
    title('N-T Plane Projection');
    legend('SV2 / Watcher', 'SV3 / Docker', 'SV1 / Chief');
    % legend('SV2 / Watcher', 'SV1 / Chief');
    
    %sgtitle(title_string, 'FontWeight', 'bold');
    % Save figure
    saveas(gcf, fig_path);
end
