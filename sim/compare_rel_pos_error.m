function compare_rel_pos_error(t, t_orbit, pos1, pos2, vel1, vel2, title_string,fig_path)
    % Compute the error (component-wise)
    error_pos = pos2 - pos1;
    error_vel = vel2 - vel1;

    % Plotting
    figure('Color', 'w');
    component_labels_pos = {'X Error [km]', 'Y Error [km]', 'Z Error [km]'};
    component_labels_vel = {'Vx Error [km/s]', 'Vy Error [km/s]', 'Vz Error [km/s]'};

    for i = 1:3
        subplot(3,2,2*i-1);
        plot(t/t_orbit, error_pos(:,i), 'r', 'LineWidth', 1.5);
        grid on;
        ylabel(component_labels_pos{i});
        if i == 3
            %xlabel('Time [s]');
            xlabel('Number of Orbits')
        end

        subplot(3, 2, 2*i);
        plot(t/t_orbit, error_vel(:,i), 'r', 'LineWidth', 1.5);
        grid on;
        ylabel(component_labels_vel{i});
        if i == 3
            %xlabel('Time [s]');
            xlabel('Number of Orbits')
        end
    end

    sgtitle(title_string, 'FontWeight', 'bold');
    saveas(gcf, fig_path);