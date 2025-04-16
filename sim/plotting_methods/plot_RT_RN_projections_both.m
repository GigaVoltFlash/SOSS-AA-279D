function plot_RT_RN_projections_both(SV2_rel_pos1, SV2_rel_pos2, SV3_rel_pos1, SV3_rel_pos2, title_string,fig_path)
    % Extract coordinates
    x1_1 = SV2_rel_pos1(:,1); y1_1 = SV2_rel_pos1(:,2); z1_1 = SV2_rel_pos1(:,3);
    x2_1 = SV3_rel_pos1(:,1); y2_1 = SV3_rel_pos1(:,2); z2_1 = SV3_rel_pos1(:,3);
    x1_2 = SV2_rel_pos2(:,1); y1_2 = SV2_rel_pos2(:,2); z1_2 = SV2_rel_pos2(:,3);
    x2_2 = SV3_rel_pos2(:,1); y2_2 = SV3_rel_pos2(:,2); z2_2 = SV3_rel_pos2(:,3);

    figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
    
    % === X-Y Plane ===
    subplot(1, 2, 1);
    plot(y1_1, x1_1, 'g', 'LineWidth', 1.5); hold on;
    plot(y1_2, x1_2, 'b--', 'LineWidth', 1.5); 
    plot(y2_1, x2_1, 'k', 'LineWidth', 1.5);
    plot(y2_2, x2_2, 'r--', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    ylabel('Radial [km]'); xlabel('Tangential [km]');
    ylim([-2.0, 2.0]);
    title('R-T Plane Projection');
    legend('SV2 / Watcher with relative EOM', 'SV2 / Watcher with FODE', 'SV3 / Docker with relative EOM', 'SV3 / Docker with FODE', 'SV1 / Chief');

    % === X-Z Plane ===
    subplot(1, 2, 2);
    plot(z1_1, x1_1, 'g', 'LineWidth', 1.5); hold on;
    plot(z1_2, x1_2, 'b--', 'LineWidth', 1.5); 
    plot(z2_1, x2_1, 'k', 'LineWidth', 1.5);
    plot(z2_2, x2_2, 'r--', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    ylabel('Radial [km]'); xlabel('Normal [km]');
    ylim([-2.0, 2.0]);
    title('R-N Plane Projection');
    legend('SV2 / Watcher with relative EOM', 'SV2 / Watcher with FODE', 'SV3 / Docker with relative EOM', 'SV3 / Docker with FODE', 'SV1 / Chief');
    
    sgtitle(title_string, 'FontWeight', 'bold');
    % Save figure
    saveas(gcf, fig_path);
end
