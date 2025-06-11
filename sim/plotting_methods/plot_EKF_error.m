function plot_EKF_error(full_times, t_orbit, error, covariance, save_path)
    num_orbits = full_times / t_orbit;
    roe_labels = {'\delta a [m]', '\delta \lambda [m]', '\delta e_x [m]', ...
                  '\delta e_y [m]', '\delta i_x [m]', '\delta i_y [m]'};

    fig = figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    h_all = [];

    for i = 1:6
        ax = nexttile(i);
        hold(ax, 'on'); grid on;

        std_1 = sqrt(covariance(:, i, i));
        std_2 = 2 * std_1;
        std_3 = 3 * std_1;
        
        % 3σ Region (orange)
        h_fill3 = fill(ax, [num_orbits; flipud(num_orbits)], ...
            [std_3; flipud(-std_3)], ...
            [1.0, 0.8, 0.6], 'FaceAlpha', 0.8, 'EdgeColor', 'none');

        % Fill ±2σ Region (light purple)
        h_fill2 = fill(ax, [num_orbits; flipud(num_orbits)], ...
            [std_2; flipud(-std_2)], [0.85 0.75 1.0], ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none');

        % Fill ±1σ Region (light red)
        h_fill1 = fill(ax, [num_orbits; flipud(num_orbits)], ...
            [std_1; flipud(-std_1)], [1.0 0.6 0.6], ...
            'FaceAlpha', 0.4, 'EdgeColor', 'none');

        % Plot EKF error
        h_err = plot(ax, num_orbits, error(:, i), 'k', 'LineWidth', 1.5);

        xlabel(ax, 'Time [orbits]');
        ylabel(ax, roe_labels{i});

        if i == 1
            h_all = [h_err, h_fill3, h_fill2, h_fill1];
        end
    end

    % Unified legend below all subplots
    lg = legend(h_all, {'EKF Error','3\sigma Region', '2\sigma Region', '1\sigma Region'}, ...
        'Orientation', 'horizontal', ...
        'Location', 'southoutside', ...
        'Box', 'off', 'FontSize', 10);

    lg.Layout.Tile = 'south';

    % Save full plot
    saveas(fig, save_path);

    % Save zoomed-in version (first orbit)
    % for i = 1:6
    %     nexttile(i);
    %     xlim([0, 1]);
    % end
    % zoomed_fig_path = strrep(save_path, '.png', '_zoomed.png');
    % saveas(fig, zoomed_fig_path);
end
