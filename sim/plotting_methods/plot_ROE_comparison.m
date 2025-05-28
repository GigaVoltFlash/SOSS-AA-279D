function plot_ROE_comparison(t, t_orbit, ROE1, ROE2, ...
    label1, label2, title_str2, fig_path2, title_str3, fig_path3)

    d_a1 = ROE1(:, 1);
    d_lambda1 = ROE1(:, 2);
    d_e_x1 = ROE1(:, 3);
    d_e_y1 = ROE1(:, 4);
    d_i_x1 = ROE1(:, 5);
    d_i_y1 = ROE1(:, 6);

    d_a2 = ROE2(:, 1);
    d_lambda2 = ROE2(:, 2);
    d_e_x2 = ROE2(:, 3);
    d_e_y2 = ROE2(:, 4);
    d_i_x2 = ROE2(:, 5);
    d_i_y2 = ROE2(:, 6);

    % === Plane Projections ===
    figure('Color', 'w', 'Position', [100, 100, 1000, 600]);

    ax1 = subplot(1,3,1); hold on; axis equal; grid on;
    xlabel('\delta e_x [m]'); ylabel('\delta e_y [m]');

    ax2 = subplot(1,3,2); hold on; axis equal; grid on;
    xlabel('\delta i_x [m]'); ylabel('\delta i_y [m]');

    ax3 = subplot(1,3,3); hold on; axis equal; grid on;
    xlabel('\delta \lambda [m]'); ylabel('\delta a [m]');

    lw = 2.0;
    plot(ax1, d_e_x1, d_e_y1, 'b', 'LineWidth', lw);
    plot(ax1, d_e_x2, d_e_y2, 'r--', 'LineWidth', lw);

    plot(ax2, d_i_x1, d_i_y1, 'b', 'LineWidth', lw);
    plot(ax2, d_i_x2, d_i_y2, 'r--', 'LineWidth', lw);

    plot(ax3, d_lambda1, d_a1, 'b', 'LineWidth', lw);
    plot(ax3, d_lambda2, d_a2, 'r--', 'LineWidth', lw);

    legend(ax1, label1, label2, 'Location', 'best');
    legend(ax2, label1, label2, 'Location', 'best');
    legend(ax3, label1, label2, 'Location', 'best');

    %sgtitle(title_str2);
    saveas(gcf, fig_path2);

    % === Time Histories ===
    time_orbits = t / t_orbit;

    figure('Color', 'w', 'Position', [100, 100, 800, 600]);

    axs = gobjects(6,1);
    labels = {'\delta a [m]', '\delta \lambda [m]', '\delta e_x [m]', ...
              '\delta e_y [m]', '\delta i_x [m]', '\delta i_y [m]'};

    ydata1 = {d_a1, d_lambda1, d_e_x1, d_e_y1, d_i_x1, d_i_y1};
    ydata2 = {d_a2, d_lambda2, d_e_x2, d_e_y2, d_i_x2, d_i_y2};

    for k = 1:6
        axs(k) = subplot(3,2,k); hold on; grid on;
        xlabel('Number of orbits'); ylabel(labels{k});
        plot(time_orbits, ydata1{k}, 'b', 'LineWidth', lw);
        plot(time_orbits, ydata2{k}, 'r--', 'LineWidth', lw);
        if k==1
            legend(label1, label2, 'Location', 'best');
        end
    end

    %sgtitle(title_str3);
    saveas(gcf, fig_path3);
end
