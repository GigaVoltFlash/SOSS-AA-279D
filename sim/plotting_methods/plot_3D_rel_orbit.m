function [] = plot_3D_rel_orbit(SV2_rel_pos1, SV3_rel_pos1, both_on_one, title_string, fig_path)
    x1_1 = SV2_rel_pos1(:,1); y1_1 = SV2_rel_pos1(:,2); z1_1 = SV2_rel_pos1(:,3);
    x2_1 = SV3_rel_pos1(:,1); y2_1 = SV3_rel_pos1(:,2); z2_1 = SV3_rel_pos1(:,3);
    
    if both_on_one == false
        figure();
    
        subplot(1, 2, 1);
        plot3(x1_1, y1_1, z1_1, 'b');
        hold on;
    
    
        xlabel('Radial [km]');
        ylabel('Tangential [km]');
        zlabel('Normal [km]');
        title('SV2 3D Relative Orbit')
        grid on;
        axis equal;
        %legend('YA Dif. Eq.', 'YA Mapping');
        hold off;
        %saveas(gcf, fig_path_SV2);
    
        subplot(1, 2, 2);
        plot3(x2_1, y2_1, z2_1, 'b');
        hold on;
    
    
        xlabel('Radial [km]');
        ylabel('Tangential [km]');
        zlabel('Normal [km]');
        title('SV3 3D Relative Orbit')
        grid on;
        axis equal;
        %legend('YA Dif. Eq.', 'YA Mapping');
        hold off;
    
        %sgtitle(title_string, 'FontWeight', 'bold');
        % Save figure
        saveas(gcf, fig_path);
    else
        figure();

        % Plot SV2 (Watcher)
        plot3(x1_1, y1_1, z1_1, 'b', 'DisplayName', 'SV2 / Watcher');
        hold on;
        
        % Plot SV3 (Docker)
        plot3(x2_1, y2_1, z2_1, 'r', 'DisplayName', 'SV3 / Docker');
        
        % Plot SV1 (Chief) as black dot at origin
        plot3(0, 0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', 'SV1 / Chief');
        
        xlabel('Radial [km]');
        ylabel('Tangential [km]');
        zlabel('Normal [km]');
        grid on;
        axis equal;
        legend('Location', 'best');
        %sgtitle(title_string, 'FontWeight', 'bold');
        
        % Save figure
        saveas(gcf, fig_path);
    end
end