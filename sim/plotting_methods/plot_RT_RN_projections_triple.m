function plot_RT_RN_projections_triple(SV2_rel_pos1, SV2_rel_pos2, SV2_rel_pos3, SV3_rel_pos1, SV3_rel_pos2, SV3_rel_pos3, title_string,fig_path)
    % Extract coordinates
    x1_1 = SV2_rel_pos1(:,1); y1_1 = SV2_rel_pos1(:,2); z1_1 = SV2_rel_pos1(:,3);
    x2_1 = SV3_rel_pos1(:,1); y2_1 = SV3_rel_pos1(:,2); z2_1 = SV3_rel_pos1(:,3);
    x1_2 = SV2_rel_pos2(:,1); y1_2 = SV2_rel_pos2(:,2); z1_2 = SV2_rel_pos2(:,3);
    x2_2 = SV3_rel_pos2(:,1); y2_2 = SV3_rel_pos2(:,2); z2_2 = SV3_rel_pos2(:,3);
    x1_3 = SV2_rel_pos3(:,1); y1_3 = SV2_rel_pos3(:,2); z1_3 = SV2_rel_pos3(:,3);
    x2_3 = SV3_rel_pos3(:,1); y2_3 = SV3_rel_pos3(:,2); z2_3 = SV3_rel_pos3(:,3);

    figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
    
    % === X-Y Plane ===
    subplot(1, 3, 1);
    plot(y1_1, x1_1, 'c', 'LineWidth', 1.5); hold on;
    plot(y1_2, x1_2, 'g', 'LineWidth', 1.5); 
    plot(y1_3, x1_3, 'b', 'LineWidth', 1.5);
    plot(y2_1, x2_1, 'm--', 'LineWidth', 1.5);
    plot(y2_2, x2_2, 'y--', 'LineWidth', 1.5);
    plot(y2_3, x2_3, 'r--', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    axis equal;
    ylabel('Radial [km]'); xlabel('Tangential [km]');
    %ylim([-2.0, 2.0]);
    title('R-T Plane Projection');
    legend('SV2 / Watcher with YA Dif. Eq.', 'SV2 / Watcher with YA Mapping', 'SV2 / Watcher with FERM',...
        'SV3 / Watcher with YA Dif. Eq.', 'SV3 / Watcher with YA Mapping', 'SV3 / Watcher with FERM',...
        'SV1 / Chief');

    % === X-Z Plane ===
    subplot(1, 3, 2);
    plot(z1_1, x1_1, 'c', 'LineWidth', 1.5); hold on;
    plot(z1_2, x1_2, 'g', 'LineWidth', 1.5); 
    plot(z1_3, x1_3, 'b', 'LineWidth', 1.5);
    plot(z2_1, x2_1, 'm--', 'LineWidth', 1.5);
    plot(z2_2, x2_2, 'y--', 'LineWidth', 1.5);
    plot(z2_3, x2_3, 'r--', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    axis equal;
    ylabel('Radial [km]'); xlabel('Normal [km]');
    %ylim([-2.0, 2.0]);
    title('R-N Plane Projection');
    legend('SV2 / Watcher with YA Dif. Eq.', 'SV2 / Watcher with YA Mapping', 'SV2 / Watcher with FERM',...
        'SV3 / Watcher with YA Dif. Eq.', 'SV3 / Watcher with YA Mapping', 'SV3 / Watcher with FERM',...
        'SV1 / Chief');    
    % === Y-Z Plane ===
    subplot(1, 3, 3);
    plot(y1_1, z1_1, 'c', 'LineWidth', 1.5); hold on;
    plot(y1_2, z1_2, 'g', 'LineWidth', 1.5); 
    plot(y1_3, z1_3, 'b', 'LineWidth', 1.5);
    plot(y2_1, z2_1, 'm--', 'LineWidth', 1.5);
    plot(y2_2, z2_2, 'y--', 'LineWidth', 1.5);
    plot(y2_3, z2_3, 'r--', 'LineWidth', 1.5);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % Origin
    grid on;
    axis equal;
    ylabel('Normal [km]'); xlabel('Tangential [km]');
    title('N-T Plane Projection');
    legend('SV2 / Watcher with YA Dif. Eq.', 'SV2 / Watcher with YA Mapping', 'SV2 / Watcher with FERM',...
        'SV3 / Watcher with YA Dif. Eq.', 'SV3 / Watcher with YA Mapping', 'SV3 / Watcher with FERM',...
        'SV1 / Chief');


    sgtitle(title_string, 'FontWeight', 'bold');
    % Save figure
    saveas(gcf, fig_path);
end
