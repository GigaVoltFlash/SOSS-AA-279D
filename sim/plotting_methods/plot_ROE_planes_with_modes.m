function plot_ROE_planes_with_modes(t, t_orbit, d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y, SV3_modes, num_orbits_modes, num_orbits_station_keep, filename1, filename2, filename3)        
    colors = lines(size(SV3_modes, 1));
    figure('Color', 'w','Position', [100, 100, 1000, 600]);
    
    ax1 = subplot(1,3,1);
    hold(ax1, 'on'); axis equal; grid on;
    xlabel('\delta e_x [m]'); ylabel('\delta e_y [m]');
    lw = 2.5;
    plot(ax1, d_e_x, d_e_y, 'k-', 'LineWidth', lw);

    ax2 = subplot(1,3,2);
    hold(ax2, 'on'); axis equal; grid on;
    xlabel('\delta i_x [m]'); ylabel('\delta i_y [m]');
    plot(ax2, d_i_x, d_i_y, 'k-', 'LineWidth', lw);

    ax3 = subplot(1,3,3);
    hold(ax3, 'on'); axis equal; grid on;
    xlabel('\delta \lambda [m]'); ylabel('\delta a [m]');
    plot(ax3, d_lambda, d_a, 'k-', 'LineWidth', lw);

    % Initialize current positions to first mode:
    current_pos_RT = SV3_modes(1,3:4);
    current_pos_RN = SV3_modes(1,5:6);
    current_pos_NR = SV3_modes(1,[2,1]);
    
    for i = 2:size(SV3_modes,1)
        color_i = colors(i,:);
        
        delta_RT = SV3_modes(i,3:4) - SV3_modes(i-1,3:4);
        delta_RN = SV3_modes(i,5:6) - SV3_modes(i-1,5:6);
        delta_NR = SV3_modes(i,[2,1]) - SV3_modes(i-1,[2,1]);
        
        % R-T plane
        quiver(ax1, current_pos_RT(1), current_pos_RT(2), delta_RT(1), delta_RT(2), 0, ...
               'Color', color_i, 'LineWidth', 1.5);
        current_pos_RT = current_pos_RT + delta_RT;
        
        % R-N plane
        quiver(ax2, current_pos_RN(1), current_pos_RN(2), delta_RN(1), delta_RN(2), 0, ...
               'Color', color_i, 'LineWidth', 1.5);
        current_pos_RN = current_pos_RN + delta_RN;
        
        % N-R plane
        quiver(ax3, current_pos_NR(1), current_pos_NR(2), delta_NR(1), delta_NR(2), 0, ...
               'Color', color_i, 'LineWidth', 1.5);
        current_pos_NR = current_pos_NR + delta_NR;
    end

    % Add unified legend below all subplots
    lgd1 = legend([arrayfun(@(i) plot(NaN,NaN,'Color',colors(i,:),'LineWidth',3), 1:size(SV3_modes,1))'], ...
        arrayfun(@(i) sprintf('Mode %d', i), 1:size(SV3_modes,1), 'UniformOutput', false));
    lgd1.Orientation = 'horizontal';
    lgd1.Units = 'normalized';
    lgd1.Position = [0.2, 0.01, 0.6, 0.03];
    lgd1.FontSize = 10;

    hold(ax1, 'off'); hold(ax2, 'off'); hold(ax3, 'off');
    saveas(gcf, filename1);

    % Time-series plot
    time_orbits = t / t_orbit;
    figure('Position', [100, 100, 800, 700]);

    ax1_t = subplot(3,2,1); hold(ax1_t, 'on'); grid on; xlabel('Orbits'); ylabel('\delta a [m]');
    ax2_t = subplot(3,2,2); hold(ax2_t, 'on'); grid on; xlabel('Orbits'); ylabel('\delta \lambda [m]');
    ax3_t = subplot(3,2,3); hold(ax3_t, 'on'); grid on; xlabel('Orbits'); ylabel('\delta e_x [m]');
    ax4_t = subplot(3,2,4); hold(ax4_t, 'on'); grid on; xlabel('Orbits'); ylabel('\delta e_y [m]');
    ax5_t = subplot(3,2,5); hold(ax5_t, 'on'); grid on; xlabel('Orbits'); ylabel('\delta i_x [m]');
    ax6_t = subplot(3,2,6); hold(ax6_t, 'on'); grid on; xlabel('Orbits'); ylabel('\delta i_y [m]');

    plot(ax1_t, time_orbits, d_a, 'k-', 'LineWidth', lw);
    plot(ax2_t, time_orbits, d_lambda, 'k-', 'LineWidth', lw);
    plot(ax3_t, time_orbits, d_e_x, 'k-', 'LineWidth', lw);
    plot(ax4_t, time_orbits, d_e_y, 'k-', 'LineWidth', lw);
    plot(ax5_t, time_orbits, d_i_x, 'k-', 'LineWidth', lw);
    plot(ax6_t, time_orbits, d_i_y, 'k-', 'LineWidth', lw);

    % Overlay desired values and shaded station-keeping
    cum_orbits = 0;
    sk_patch = [0.9 0.9 0.9];
    mode_handles = gobjects(size(SV3_modes, 1), 1);
    for i = 1:size(SV3_modes, 1)
        if i > length(num_orbits_modes) || i > length(num_orbits_station_keep)
            break;
        end
        orb_start = cum_orbits;
        orb_end = orb_start + num_orbits_modes(i);
        cum_orbits = orb_end + num_orbits_station_keep(i);

        color_i = colors(i,:);
        mode_handles(i) = plot(ax1_t, [orb_start orb_end], SV3_modes(i,1)*[1 1], '-', 'Color', color_i, 'LineWidth', 2);
        plot(ax2_t, [orb_start orb_end], SV3_modes(i,2)*[1 1], '-', 'Color', color_i, 'LineWidth', 2);
        plot(ax3_t, [orb_start orb_end], SV3_modes(i,3)*[1 1], '-', 'Color', color_i, 'LineWidth', 2);
        plot(ax4_t, [orb_start orb_end], SV3_modes(i,4)*[1 1], '-', 'Color', color_i, 'LineWidth', 2);
        plot(ax5_t, [orb_start orb_end], SV3_modes(i,5)*[1 1], '-', 'Color', color_i, 'LineWidth', 2);
        plot(ax6_t, [orb_start orb_end], SV3_modes(i,6)*[1 1], '-', 'Color', color_i, 'LineWidth', 2);

        % Add shaded SK patch
        sk_start = orb_end;
        sk_end = cum_orbits;
        for ax = [ax1_t, ax2_t, ax3_t, ax4_t, ax5_t, ax6_t]
            ylimits = get(ax, 'YLim');
            fill(ax, [sk_start sk_end sk_end sk_start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], sk_patch, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        end
    end

    % Add unified legend using handles, positioned below all subplots
    lgd = legend([mode_handles(1:end); plot(NaN,NaN,'Color',sk_patch)], ...
           [arrayfun(@(i) sprintf('Mode %d', i), 1:size(SV3_modes,1), 'UniformOutput', false), {'Station-Keeping'}]);
    lgd.Orientation = 'horizontal';
    lgd.Units = 'normalized';
    lgd.Position = [0.1, 0.01, 0.8, 0.03];
    lgd.FontSize = 10;

    hold(ax1_t, 'off'); hold(ax2_t, 'off'); hold(ax3_t, 'off');
    hold(ax4_t, 'off'); hold(ax5_t, 'off'); hold(ax6_t, 'off');
    saveas(gcf, filename2);

    % === ROE ERROR PLOTS ===
    d_err = [d_a(:), d_lambda(:), d_e_x(:), d_e_y(:), d_i_x(:), d_i_y(:)];
    err = zeros(size(d_err));

    cum_orbits = 0;
    time_orbits = t / t_orbit;
    prev_desired = zeros(1,6);

    for i = 1:length(SV3_modes)
        if i > length(num_orbits_modes) || i > length(num_orbits_station_keep)
            break;
        end
        idx_start = find(time_orbits >= cum_orbits, 1, 'first');
        idx_end = find(time_orbits >= cum_orbits + num_orbits_modes(i), 1, 'first');
        if isempty(idx_end), idx_end = length(t); end

        for j = 1:6
            err(idx_start:idx_end, j) = abs(d_err(idx_start:idx_end, j) - SV3_modes(i, j));
        end

        % SK error relative to prior desired value
        sk_start = idx_end;
        sk_end = find(time_orbits >= cum_orbits + num_orbits_modes(i) + num_orbits_station_keep(i), 1, 'first');
        if isempty(sk_end), sk_end = length(t); end

        for j = 1:6
            err(sk_start:sk_end, j) = abs(d_err(sk_start:sk_end, j) - SV3_modes(i, j));
        end

        cum_orbits = cum_orbits + num_orbits_modes(i) + num_orbits_station_keep(i);
    end

    figure('Position', [100, 100, 800, 700]);
    labels = {'\delta a error [m]', '\delta \lambda error [m]', '\delta e_x error [m]', '\delta e_y error [m]', '\delta i_x error [m]', '\delta i_y error [m]'};

    for j = 1:6
        ax_err = subplot(3,2,j);
        hold(ax_err, 'on'); grid on;
        xlabel('Orbits'); ylabel(labels{j});
        plot(ax_err, time_orbits, err(:,j), 'k-', 'LineWidth', 1.5);

        cum_orbits = 0;
        for i = 1:size(SV3_modes,1)
            if i > length(num_orbits_modes) || i > length(num_orbits_station_keep)
                break;
            end
            orb_start = cum_orbits;
            orb_end = orb_start + num_orbits_modes(i);
            cum_orbits = orb_end + num_orbits_station_keep(i);

            ylimits = get(ax_err, 'YLim');
            fill(ax_err, [orb_start orb_end orb_end orb_start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], colors(i,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
            fill(ax_err, [orb_end cum_orbits cum_orbits orb_end], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], sk_patch, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        end
    end

    lgd3 = legend([arrayfun(@(i) plot(NaN,NaN,'Color',colors(i,:),'LineWidth',3), 1:size(SV3_modes,1))'; plot(NaN,NaN,'Color',sk_patch)], ...
           [arrayfun(@(i) sprintf('Mode %d', i), 1:size(SV3_modes,1), 'UniformOutput', false), {'Station-Keeping'}]);
    lgd3.Orientation = 'horizontal';
    lgd3.Units = 'normalized';
    lgd3.Position = [0.2, 0.01, 0.6, 0.03];
    lgd3.FontSize = 10;

    saveas(gcf, filename3);
end
