function compare_rel_pos_error(t, pos1, pos2, vel1, vel2, fig_path)
    % Compute the error (component-wise)
    error_pos = pos2 - pos1;
    error_vel = vel2 - vel1;

    % Plotting
    figure('Color', 'w');
    component_labels_pos = {'X Error [km]', 'Y Error [km]', 'Z Error [km]'};
    component_labels_vel = {'Vx Error [km/s]', 'Vy Error [km/s]', 'Vz Error [km/s]'};

    for i = 1:3
        subplot(3,2,2*i-1);
        plot(t, error_pos(:,i), 'r', 'LineWidth', 1.5);
        grid on;
        ylabel(component_labels_pos{i});
        if i == 3
            xlabel('Time [s]');
        end

        subplot(3, 2, 2*i);
        plot(t, error_vel(:,i), 'r', 'LineWidth', 1.5);
        grid on;
        ylabel(component_labels_vel{i});
        if i == 3
            xlabel('Time [s]');
        end
    end

    sgtitle('Relative Position & Velocity Error Components Over Time');
    saveas(gcf, fig_path);