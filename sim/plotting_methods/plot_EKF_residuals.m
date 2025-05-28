function plot_EKF_residuals(full_times, t_orbit, pre_fit, post_fit, noise, save_path)
    num_orbits = full_times / t_orbit;
    roe_labels = {'\delta a [m]', '\delta \lambda [m]', '\delta e_x [m]', ...
                  '\delta e_y [m]', '\delta i_x [m]', '\delta i_y [m]'};
    
    fig = figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % For legend handles (create only once)
    h_all = [];

    for i = 1:6
        ax = nexttile(i);
        hold(ax, 'on');
        grid(ax, 'on');

        % Scatter plots
        h_noise = scatter(ax, num_orbits, noise(:,i), 10, 'k', 'filled', ...
            'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.3);
        h_pre   = scatter(ax, num_orbits, pre_fit(:,i), 10, 'b', 'filled', ...
            'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', 0.7);
        h_post  = scatter(ax, num_orbits, post_fit(:,i), 8, 'r', '^', 'filled', ...
            'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', 1);

        xlabel(ax, 'Time [orbits]');
        ylabel(ax, roe_labels{i});
        xlim(ax, [0, 1]);

        % Store legend handles only on first tile
        if i == 1
            h_all = [h_noise, h_pre, h_post];
        end
    end

    % Add unified legend below all subplots
    lg = legend(h_all, {'Measurement Noise', 'Pre-fit Residuals', 'Post-fit Residuals'}, ...
        'Orientation', 'horizontal', ...
        'Location', 'southoutside', ...
        'Box', 'off', 'FontSize', 10);

    lg.Layout.Tile = 'south';

    saveas(fig, save_path);
end
