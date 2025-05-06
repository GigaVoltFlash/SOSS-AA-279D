function plot_RT_RN_projections_separate(SV2_rel_pos, SV3_rel_pos, RT, RN, NT, title_string, fig_path_base)
    % Extract coordinates (km --> m)
    x1 = SV2_rel_pos(:,1)*1e3; y1 = SV2_rel_pos(:,2)*1e3; z1 = SV2_rel_pos(:,3)*1e3;
    x2 = SV3_rel_pos(:,1)*1e3; y2 = SV3_rel_pos(:,2)*1e3; z2 = SV3_rel_pos(:,3)*1e3;

    % === R-T Plane ===
    if RT
        figure('Color', 'w');
        plot(y1, x1, 'b', 'LineWidth', 1.5); hold on;
        plot(y2, x2, 'r', 'LineWidth', 1.5);
        plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
        grid on;
        axis equal;
        ylabel('Radial [m]'); xlabel('Tangential [m]');
        %title('R-T Plane Projection');
        legend('SV2 / Watcher', 'SV3 / Docker', 'SV1 / Chief');
        sgtitle([title_string ' — R-T Plane']);
        saveas(gcf, [fig_path_base '_RT.png']);
    end

    % === R-N Plane ===
    if RN
        figure('Color', 'w');
        plot(z1, x1, 'b', 'LineWidth', 1.5); hold on;
        plot(z2, x2, 'r', 'LineWidth', 1.5);
        plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
        grid on;
        axis equal;
        ylabel('Radial [m]'); xlabel('Normal [m]');
        %title('R-N Plane Projection');
        legend('SV2 / Watcher', 'SV3 / Docker', 'SV1 / Chief');
        sgtitle([title_string ' — R-N Plane']);
        saveas(gcf, [fig_path_base '_RN.png']);
    end

    % === N-T Plane ===
    if NT
        figure('Color', 'w');
        plot(y1, z1, 'b', 'LineWidth', 1.5); hold on;
        plot(y2, z2, 'r', 'LineWidth', 1.5);
        plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
        grid on;
        axis equal;
        ylabel('Normal [m]'); xlabel('Tangential [m]');
        %title('N-T Plane Projection');
        legend('SV2 / Watcher', 'SV3 / Docker', 'SV1 / Chief');
        sgtitle([title_string ' — N-T Plane']);
        saveas(gcf, [fig_path_base '_NT.png']);
    end
end
