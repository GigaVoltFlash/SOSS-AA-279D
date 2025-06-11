function plot_delta_v_cumulative_timeline(delta_v_times, delta_v_vals, t_orbit, SV3_modes, num_orbits_modes, num_orbits_station_keep, is_continuous, filename)
    % Converts delta-v values (in RTN m/s) into a time series plot of R, T, N components with
    % shaded regions for operational modes and station-keeping, and plots cumulative delta-v

    colors = lines(size(SV3_modes, 1));
    time_orbits = delta_v_times / t_orbit;

    figure('Position', [100, 100, 800, 700]);
    ax = gobjects(4,1);
    labels = {'Radial [m/s]', 'Tangential [m/s]', 'Normal [m/s]', 'Cumulative \Deltav [m/s]'};

    for i = 1:4
        ax(i) = subplot(4,1,i);
        hold(ax(i), 'on'); grid on;
        ylabel(ax(i), labels{i});
    end
    xlabel(ax(4), 'Time [orbits]');

    % Plot impulses or continuous values
    if is_continuous
        plot(ax(1), time_orbits, delta_v_vals(:,1));
        plot(ax(2), time_orbits, delta_v_vals(:,2));
        plot(ax(3), time_orbits, delta_v_vals(:,3));
    else
        stem(ax(1), time_orbits, delta_v_vals(:,1), 'filled', 'LineWidth', 1.5);
        stem(ax(2), time_orbits, delta_v_vals(:,2), 'filled', 'LineWidth', 1.5);
        stem(ax(3), time_orbits, delta_v_vals(:,3), 'filled', 'LineWidth', 1.5);
    end

    % Compute cumulative delta-v magnitude and plot
    dv_abs_sum = sum(abs(delta_v_vals), 2);  % sum of absolute value of each entry per row
    dv_cum = cumsum(dv_abs_sum);             % cumulative sum of that
    
    plot(ax(4), time_orbits, dv_cum, 'k', 'LineWidth', 1.5);

    % Add shaded regions and legend entries
    cum_orbits = 0;
    sk_patch = [0.9 0.9 0.9];
    mode_handles = gobjects(size(SV3_modes,1), 1);
    total_dv_per_mode = zeros(size(SV3_modes,1), 1);

    for i = 1:size(SV3_modes, 1)
        orb_start = cum_orbits;
        orb_end = orb_start + num_orbits_modes(i);
        cum_orbits = orb_end + num_orbits_station_keep(i);

        color_i = colors(i,:);
        for j = 1:4
            ylimits = get(ax(j), 'YLim');
            fill(ax(j), [orb_start orb_end orb_end orb_start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], color_i, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
            fill(ax(j), [orb_end cum_orbits cum_orbits orb_end], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], sk_patch, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        end

        % Compute total delta-v in this mode
        t_start = orb_start * t_orbit;
        t_end = orb_end * t_orbit;
        idx = delta_v_times >= t_start & delta_v_times < t_end;
        total_dv_per_mode(i) = sum(dv_abs_sum(idx));
    end

    % Print total delta-v per mode
    fprintf('\n=== Total Δv per mode ===\n');
    for i = 1:length(total_dv_per_mode)
        fprintf('Mode %d: %.6f m/s\n', i, total_dv_per_mode(i)); 
    end
    fprintf('==========================\n\n');

    % Legend
    lgd = legend([arrayfun(@(i) plot(NaN,NaN,'Color',colors(i,:),'LineWidth',3), 1:size(SV3_modes,1))'; plot(NaN,NaN,'Color',sk_patch)], ...
           [arrayfun(@(i) sprintf('Mode %d', i), 1:size(SV3_modes,1), 'UniformOutput', false), {'Station-Keeping'}]);
    lgd.Orientation = 'horizontal';
    lgd.Units = 'normalized';
    lgd.Position = [0.2, 0.01, 0.6, 0.03];
    lgd.FontSize = 10;

    saveas(gcf, filename);
end
