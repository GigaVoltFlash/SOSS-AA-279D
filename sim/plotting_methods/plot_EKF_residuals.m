function plot_EKF_residuals(full_times, t_orbit, pre_fit, post_fit, noise, save_path)
    num_orbits = full_times / t_orbit;
    figure;
    hold on;
    roe_labels = {'\delta a/a', '\delta \lambda', '\delta e_x', '\delta e_y', '\delta i_x', '\delta i_y'};
    for i = 1:6
        subplot(3,2,i);
        hold on;
        scatter(num_orbits, noise(:,i), 20, 'k', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Measurement Noise');
        scatter(num_orbits, pre_fit(:,i), 36, 'b', 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Pre-fit Residuals');
        scatter(num_orbits, post_fit(:,i), 50, 'r', '^', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', 1, 'DisplayName', 'Post-fit Residuals');
        % scatter(num_orbits, noise(:,i), 'black', 'DisplayName', 'Measurement Noise');
        % scatter(num_orbits, pre_fit(:,i), 'blue', 'DisplayName', 'Pre-fit Residuals');
        % scatter(num_orbits, post_fit(:,i), 'red', 'DisplayName', 'Post-fit Residuals');
        xlim([0, 5]);
        xlabel('Time [orbits]');
        ylabel(roe_labels(i));
        if i == 1
            legend('show', 'Location', 'best');
        end
        grid on;
        hold off;
    end
    saveas(gcf, save_path);
end