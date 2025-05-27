function plot_EKF_error(full_times, t_orbit, error,save_path)
    num_orbits = full_times / t_orbit;
    figure;
    hold on;
    roe_labels = {'\delta a/a', '\delta \lambda', '\delta e_x', '\delta e_y', '\delta i_x', '\delta i_y'};
    for i = 1:6
        subplot(3,2,i);
        hold on;
        plot(num_orbits, error(:,i), 'k:', 'LineWidth', 1.5, 'DisplayName', 'EKF Error');
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