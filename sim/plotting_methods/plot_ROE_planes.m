function plot_ROE_planes(t, t_orbit, d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y, filename1, filename2)        
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
    
    saveas(gcf, filename1);

    time_orbits = t / t_orbit;
    % 3. Relative Elements (Mean)
    figure('Position', [100, 100, 800, 600]);
     
    ax1_t = subplot(3,2,1);
    hold(ax1_t, 'on');
    grid on;
    xlabel('Number of orbits');
    ylabel('\delta a [m]');

    ax2_t = subplot(3,2,2);
    hold(ax2_t, 'on');
    grid on;
    xlabel('Number of orbits');
    ylabel('\delta \lambda [m]');

    ax3_t = subplot(3,2,3);
    hold(ax3_t, 'on');
    grid on;
    xlabel('Number of orbits');
    ylabel('\delta e_x [m]');

    ax4_t = subplot(3,2,4);
    hold(ax4_t, 'on');
    grid on;
    xlabel('Number of orbits');
    ylabel('\delta e_y [m]');

    ax5_t = subplot(3,2,5);
    hold(ax5_t, 'on');
    grid on;
    xlabel('Number of orbits');
    ylabel('\delta i_x [m]');

    ax6_t = subplot(3,2,6);
    hold(ax6_t, 'on');
    grid on;
    xlabel('Number of orbits');
    ylabel('\delta i_y [m]');
    
    lw = 2.5;
    
    plot(ax1_t, time_orbits, d_a, 'LineWidth', lw);
    plot(ax2_t, time_orbits, d_lambda, 'LineWidth', lw);
    plot(ax3_t, time_orbits, d_e_x, 'LineWidth', lw);
    plot(ax4_t, time_orbits, d_e_y, 'LineWidth', lw);
    plot(ax5_t, time_orbits, d_i_x, 'LineWidth', lw);
    plot(ax6_t, time_orbits, d_i_y, 'LineWidth', lw);
    
    hold(ax1_t, 'off');
    hold(ax2_t, 'off');
    hold(ax3_t, 'off');
    hold(ax4_t, 'off');
    hold(ax5_t, 'off');
    hold(ax6_t, 'off');
    
    % lgd = legend('Orientation', 'horizontal');
    % lgd.Units = 'normalized';
    % lgd.Position = [0.35, 0.02, 0.3, 0.03]; 
    % lgd.FontSize = 10;
    
    saveas(gcf, filename2);

end