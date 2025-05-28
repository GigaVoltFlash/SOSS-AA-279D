function plot_EKF_error(full_times, t_orbit, error, covariance, save_path)
    num_orbits = full_times / t_orbit;
    figure;
    hold on;
    roe_labels = {'\delta a [m]', '\delta \lambda [m]', '\delta e_x [m]', '\delta e_y [m]', '\delta i_x [m]', '\delta i_y [m]'};
    for i = 1:6
        subplot(3,2,i);
        hold on;
        plot(num_orbits, error(:,i), 'k', 'LineWidth', 1.5, 'DisplayName', 'EKF Error');
        % Plot 1 sigma covariance bounds
        sigma = sqrt(covariance(:,i, i));
        fill([num_orbits; flipud(num_orbits)], [sigma; flipud(-sigma)], [0.7 0.85 1], ...
            'EdgeColor', [0.2 0.4 0.8], 'FaceAlpha', 0.5, 'DisplayName', '1\sigma Region');
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