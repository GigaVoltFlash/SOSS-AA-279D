function plot_EKF_residuals(full_times, t_orbit, pre_fit, post_fit, noise, save_path)
    num_orbits = full_times / t_orbit;
    % New labels: RTN xyz (top row), ECI xyz (bottom row)
    new_labels = {'RTN x', 'RTN y', 'RTN z', 'ECI x', 'ECI y', 'ECI z'};
    
    fig = figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % For legend handles (create only once)
    h_all = [];

    % Indices for RTN and ECI
    idx_order = [1 2 3 4 5 6]; % RTN: 1-3, ECI: 4-6
    for i = 1:6
        ax = nexttile(i);
        hold(ax, 'on');
        grid(ax, 'on');

        % Scatter plots
        h_noise = scatter(ax, num_orbits, noise(:,idx_order(i)), 10, 'k', 'filled', ...
            'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.3);
        h_pre   = scatter(ax, num_orbits, pre_fit(:,idx_order(i)), 10, 'b', 'filled', ...
            'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', 0.7);
        h_post  = scatter(ax, num_orbits, post_fit(:,idx_order(i)), 8, 'r', '^', 'filled', ...
            'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', 1);

        xlabel(ax, 'Time [orbits]');
        ylabel(ax, new_labels{i});
        %xlim(ax, [0, 1]);

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
