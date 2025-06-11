function plot_ROE_comparison_with_cov(t, t_orbit, ROE1, ROE2, covariance, ...
    label1, label2, title_str2, fig_path2, title_str3, fig_path3)

    % === Preprocess ===
    d_a1 = ROE1(:, 1); d_lambda1 = ROE1(:, 2);
    d_e_x1 = ROE1(:, 3); d_e_y1 = ROE1(:, 4);
    d_i_x1 = ROE1(:, 5); d_i_y1 = ROE1(:, 6);

    d_a2 = ROE2(:, 1); d_lambda2 = ROE2(:, 2);
    d_e_x2 = ROE2(:, 3); d_e_y2 = ROE2(:, 4);
    d_i_x2 = ROE2(:, 5); d_i_y2 = ROE2(:, 6);

    ydata1 = {d_a1, d_lambda1, d_e_x1, d_e_y1, d_i_x1, d_i_y1};
    ydata2 = {d_a2, d_lambda2, d_e_x2, d_e_y2, d_i_x2, d_i_y2};
    labels = {'\delta a [m]', '\delta \lambda [m]', '\delta e_x [m]', ...
              '\delta e_y [m]', '\delta i_x [m]', '\delta i_y [m]'};

    time_orbits = t / t_orbit;
    lw = 2.0;

    % === Tiled Layout ===
    fig = figure('Color', 'w', 'Position', [100, 100, 900, 700]);
    tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    h_all = []; % Store handles for unified legend

    for k = 1:6
        ax = nexttile(k);
        hold(ax, 'on'); grid on;

        std_1 = sqrt(covariance(:, k, k));
        std_2 = 2 * std_1;
        std_3 = 3 * std_1;
        
        % 3σ Region (orange)
        h_fill3 = fill(ax, [time_orbits; flipud(time_orbits)], ...
            [ROE2(:, k) + std_3; flipud(ROE2(:, k) - std_3)], ...
            [1.0, 0.8, 0.6], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        
        % 2σ Region (purple)
        h_fill2 = fill(ax, [time_orbits; flipud(time_orbits)], ...
            [ROE2(:, k) + std_2; flipud(ROE2(:, k) - std_2)], ...
            [0.85 0.75 1.0], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        
        % 1σ Region (red)
        h_fill1 = fill(ax, [time_orbits; flipud(time_orbits)], ...
            [ROE2(:, k) + std_1; flipud(ROE2(:, k) - std_1)], ...
            [1.0 0.6 0.6], 'FaceAlpha', 0.8, 'EdgeColor', 'none');

        % Plot lines
        h1 = plot(ax, time_orbits, ydata1{k}, 'b', 'LineWidth', lw);
        h2 = plot(ax, time_orbits, ydata2{k}, 'r--', 'LineWidth', lw);

        xlabel(ax, 'Number of orbits');
        ylabel(ax, labels{k});

        if k == 1
            h_all = [h1, h2, h_fill3, h_fill2, h_fill1];
        end
    end

    lg = legend(h_all, ...
        {label1, label2,'3\sigma Region', '2\sigma Region', '1\sigma Region'}, ...
        'Orientation', 'horizontal', ...
        'Location', 'southoutside', ...
        'Box', 'off', 'FontSize', 10);

    lg.Layout.Tile = 'south';

    % === Save ===
    saveas(fig, fig_path3);

    % === Save Zoomed-In Version ===
    % for k = 1:6
    %     nexttile(k);
    %     xlim([0, 1]);
    % end
    % zoomed_fig_path = strrep(fig_path3, '.png', '_zoomed.png');
    % saveas(fig, zoomed_fig_path);
end
