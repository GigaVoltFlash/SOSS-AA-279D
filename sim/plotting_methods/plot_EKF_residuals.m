function plot_EKF_residuals(full_times, t_orbit, pre_fit, post_fit, noise, save_path)
    num_orbits = full_times / t_orbit;
    figure;
    hold on;
    roe_labels = {'\delta a/a', '\delta \lambda', '\delta e_x', '\delta e_y', '\delta i_x', '\delta i_y'};
    for i = 1:6
        subplot(3,2,i);
        hold on;
        plot(num_orbits, pre_fit(:,i), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Pre-fit Residuals');
        plot(num_orbits, post_fit(:,i), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Post-fit Residuals');
        plot(num_orbits, noise(:,i), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Measurement Noise');
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