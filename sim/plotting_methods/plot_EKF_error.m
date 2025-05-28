function plot_EKF_error(full_times, t_orbit, error, covariance, save_path)
    num_orbits = full_times / t_orbit;
    figure;
    hold on;
    roe_labels = {'\delta a [m]', '\delta \lambda [m]', '\delta e_x [m]', '\delta e_y [m]', '\delta i_x [m]', '\delta i_y [m]'};
    for i = 1:6
        subplot(3,2,i);
        hold on;
        plot(num_orbits, error(:,i), 'k', 'LineWidth', 1.5, 'DisplayName', 'EKF Error');
        % Plot sigma covariance bounds

        % std_1 = sqrt(covariance(:, k, k));
        % std_2 = 2 * std_1;
        % 
        % % 2σ Region (purple)
        % fill([time_orbits; flipud(time_orbits)], ...
        %      [ROE2(:, k) + std_2; flipud(ROE2(:, k) - std_2)], ...
        %      [0.85 0.75 1.0], ...      % light purple fill 
        %      'FaceAlpha', 0.3, ...
        %      'EdgeColor', 'none');
        % 
        % % 1σ Region (red)
        % fill([time_orbits; flipud(time_orbits)], ...
        %      [ROE2(:, k) + std_1; flipud(ROE2(:, k) - std_1)], ...
        %      [1.0 0.6 0.6],'FaceAlpha', 0.4, ...
        %      'EdgeColor', 'none');
        std_1 = sqrt(covariance(:, i, i));
        std_2 = 2 * std_1;
        fill([num_orbits; flipud(num_orbits)], [std_2; flipud(-std_2)],[0.85 0.75 1.0], ...    
             'FaceAlpha', 0.3,'EdgeColor', 'none', 'DisplayName', '2\sigma Region');

        fill([num_orbits; flipud(num_orbits)], [std_1; flipud(-std_1)], [1.0 0.6 0.6],'FaceAlpha', 0.4, ...
             'EdgeColor', 'none', 'DisplayName', '1\sigma Region');
        xlabel('Time [orbits]');
        ylabel(roe_labels(i));
        if i == 1
            legend('show', 'Location', 'best');
        end
        grid on;
        hold off;
    end
    saveas(gcf, save_path);
    % Save a zoomed-in version focused on the first orbit (for convergence)
    for k = 1:6
        subplot(3,2,k);
        xlim([0, 1]);
    end
    zoomed_fig_path = strrep(save_path, '.png', '_zoomed.png');
    saveas(gcf, zoomed_fig_path);
end