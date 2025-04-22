function plot_rel_pos_vel_single(t, t_orbit, rel_pos1, rel_vel1, title_string, fig_path)
    figure('Color', 'w');
    
    % --- Relative Position Plots ---
    subplot(3, 2, 1);
    plot(t/t_orbit, rel_pos1(:,1), 'b'); grid on; hold on;
    ylabel('x [km]');
    title('Relative Position');
    
    subplot(3, 2, 3);
    plot(t/t_orbit, rel_pos1(:,2), 'b'); grid on; hold on;
    ylabel('y [km]');
    
    subplot(3, 2, 5);
    plot(t/t_orbit, rel_pos1(:,3), 'b'); grid on; hold on;
    ylabel('z [km]');
    xlabel('Number of Orbits');

    % --- Relative Velocity Plots ---
    subplot(3, 2, 2);
    plot(t/t_orbit, rel_vel1(:,1), 'b'); grid on; hold on;
    ylabel('vx [km/s]');
    title('Relative Velocity');

    subplot(3, 2, 4);
    plot(t/t_orbit, rel_vel1(:,2), 'b'); grid on; hold on;
    ylabel('vy [km/s]');

    subplot(3, 2, 6);
    plot(t/t_orbit, rel_vel1(:,3), 'b'); grid on; hold on;
    ylabel('vz [km/s]');
    %xlabel('Time [s]');
    xlabel('Number of Orbits');

    sgtitle(title_string, 'FontWeight', 'bold');
    saveas(gcf, fig_path);
end
