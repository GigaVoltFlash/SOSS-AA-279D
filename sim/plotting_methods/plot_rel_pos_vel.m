function plot_rel_pos_vel(t, t_orbit, rel_pos, rel_vel, title_string, fig_path)
    figure('Color', 'w');
    
    % --- Relative Position Plots ---
    subplot(3, 2, 1);
    plot(t/t_orbit, rel_pos(:,1), 'b'); grid on;
    ylabel('x [km]');
    title('Relative Position');
    
    subplot(3, 2, 3);
    plot(t/t_orbit, rel_pos(:,2), 'b'); grid on;
    ylabel('y [km]');
    
    subplot(3, 2, 5);
    plot(t/t_orbit, rel_pos(:,3), 'b'); grid on;
    ylabel('z [km]');
    %xlabel('Time [s]');
    xlabel('Number of Orbits');

    % --- Relative Velocity Plots ---
    subplot(3, 2, 2);
    plot(t/t_orbit, rel_vel(:,1), 'r'); grid on;
    ylabel('vx [km/s]');
    title('Relative Velocity');

    subplot(3, 2, 4);
    plot(t/t_orbit, rel_vel(:,2), 'r'); grid on;
    ylabel('vy [km/s]');

    subplot(3, 2, 6);
    plot(t/t_orbit, rel_vel(:,3), 'r'); grid on;
    ylabel('vz [km/s]');
    %xlabel('Time [s]');
    xlabel('Number of Orbits');

    sgtitle(title_string, 'FontWeight', 'bold');
    saveas(gcf, fig_path);
end
