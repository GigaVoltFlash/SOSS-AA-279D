function plot_ROE_planes(d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y, filename)        
    % Mean
    figure('Color', 'w','Position', [100, 100, 1000, 600]);
    %sgtitle(sprintf('Relative Orbital Elements (Mean) - %s', deputy_type));
    
    ax1 = subplot(1,3,1);
    hold(ax1, 'on');
    axis equal;
    grid on;
    xlabel('\delta e_x [m]');
    ylabel('\delta e_y [m]');
    %title('Relative Eccentricity Vector (Mean)');
    
    ax2 = subplot(1,3,2);
    hold(ax2, 'on');
    axis equal;
    grid on;
    xlabel('\delta i_x [m]');
    ylabel('\delta i_y [m]');
    %title('Relative Inclination Vector (Mean)');
    
    ax3 = subplot(1,3,3);
    hold(ax3, 'on');
    axis equal;
    grid on;
    xlabel('\delta \lambda [m]');
    ylabel('\delta a [m]');
    %title('Relative Mean Longitude vs Semi-Major Axis (Mean)');
    
    
    lw = 2.5;
    
    plot(ax1, d_e_x, d_e_y, 'LineWidth', lw);
    plot(ax2, d_i_x, d_i_y, 'LineWidth', lw);
    plot(ax3, d_lambda, d_a, 'LineWidth', lw);
    
    hold(ax1, 'off');
    hold(ax2, 'off');
    hold(ax3, 'off');
    
    lgd = legend('Orientation', 'horizontal');
    lgd.Units = 'normalized';
    lgd.Position = [0.35, 0.02, 0.3, 0.03]; 
    lgd.FontSize = 10;
    
    saveas(gcf, filename);
end