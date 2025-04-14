function compare_rel_pos_error(t, pos1, pos2)
    % Inputs:
    % t     - nx1 time vector
    % pos1  - nx3 array of [x y z] positions (e.g., reference or truth)
    % pos2  - nx3 array of [x y z] positions (e.g., estimated or perturbed)

    % Compute the error (component-wise)
    error = pos2 - pos1;

    % Plotting
    figure('Color', 'w');
    component_labels = {'X Error [km]', 'Y Error [km]', 'Z Error [km]'};

    for i = 1:3
        subplot(3,1,i);
        plot(t, error(:,i), 'r', 'LineWidth', 1.5);
        grid on;
        ylabel(component_labels{i});
        if i == 3
            xlabel('Time [s]');
        end
    end

    sgtitle('Relative Position Error Components Over Time');
end