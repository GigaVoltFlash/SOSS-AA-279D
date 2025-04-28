function plot_OE_ROE_mean_osc(roe_results, tstart, tint, tend, t_orbit)
% plot_OE_ROE_mean_osc
% Inputs:
%   roe_results : struct containing ROE data (organized with SV2 and SV3)
%   tstart, tint, tend : simulation start time, timestep, and end time [s]

    % Build time vector
    time = tstart:tint:tend;
    time_orbits = time / t_orbit;
    
    % Cases and J2 conditions
    cases_to_plot = {'initial2'}; % , 'initial1'
    j2_conditions = {'no_j2', 'with_j2'};

    % Legend entries you want
    custom_legends = {'IC 2 without J2', 'IC 2 with J2'}; % 'IC 1 without J2', 'IC 1 with J2',

    % Deputies
    all_deputies = {'SV2', 'SV3'};

    % Mapping of initial case to actual deputy names
    deputy_mapping = struct( ...
        'initial1', struct('SV2', 'SV2_1', 'SV3', 'SV3_1'), ...
        'initial2', struct('SV2', 'SV2_2', 'SV3', 'SV3_2') ...
    );

    % Field names and corresponding y-axis labels
    abs_fields = {'a', 'e_x', 'e_y', 'i', 'RAAN', 'u'};
    abs_labels = {'a [km]', 'e_x', 'e_y', 'i [deg]', 'RAAN [deg]', 'u [deg]'};

    rel_fields_osc = {'d_a_osc', 'd_lambda_osc', 'd_e_x_osc', 'd_e_y_osc', 'd_i_x_osc', 'd_i_y_osc'};
    rel_labels_osc = {'\Delta a [m]', '\Delta \lambda [m]', '\Delta e_x [m]', '\Delta e_y [m]', '\Delta i_x [m]', '\Delta i_y [m]'};

    rel_fields_mean = {'d_a_mean', 'd_lambda_mean', 'd_e_x_mean', 'd_e_y_mean', 'd_i_x_mean', 'd_i_y_mean'};
    rel_labels_mean = {'\Delta a [m]', '\Delta \lambda [m]', '\Delta e_x [m]', '\Delta e_y [m]', '\Delta i_x [m]', '\Delta i_y [m]'};

    % === Now for each Deputy separately (SV2, SV3) ===
    for deputy_idx = 1:length(all_deputies)
        deputy_type = all_deputies{deputy_idx};

        % 1. Absolute Elements (Mean)
        figure;
        sgtitle(sprintf('Absolute Orbital Elements (Mean) - %s', deputy_type));
        for idx = 1:length(abs_fields)
            field_name = abs_fields{idx};
            subplot(3,2,idx);
            hold on;
            curve_idx = 0; % For matching legend entries
            for init_idx = 1:length(cases_to_plot)
                for j2_idx = 1:length(j2_conditions)
                    curve_idx = curve_idx + 1;
                    init_case = cases_to_plot{init_idx};
                    j2_case = j2_conditions{j2_idx};
                    deputy_name = deputy_mapping.(init_case).(deputy_type);

                    values = roe_results.(init_case).(j2_case).(deputy_name).([field_name, '_mean']);

                    plot(time_orbits, values, 'DisplayName', custom_legends{curve_idx});
                end
            end
            %title(strrep(field_name, '_', '\_'));
            xlabel('Number of Orbits');
            ylabel(abs_labels{idx});
            grid on;
            hold off;
        end
        legend('Location', 'bestoutside');

        % 2. Absolute Elements (Oscillating)
        figure;
        sgtitle(sprintf('Absolute Orbital Elements (Oscullating) - %s', deputy_type));
        for idx = 1:length(abs_fields)
            field_name = abs_fields{idx};
            subplot(3,2,idx);
            hold on;
            curve_idx = 0;
            for init_idx = 1:length(cases_to_plot)
                for j2_idx = 1:length(j2_conditions)
                    curve_idx = curve_idx + 1;
                    init_case = cases_to_plot{init_idx};
                    j2_case = j2_conditions{j2_idx};
                    deputy_name = deputy_mapping.(init_case).(deputy_type);

                    values = roe_results.(init_case).(j2_case).(deputy_name).([field_name, '_osc']);

                    plot(time_orbits, values, 'DisplayName', custom_legends{curve_idx});
                end
            end
            %title(strrep(field_name, '_', '\_'));
            xlabel('Number of Orbits');
            ylabel(abs_labels{idx});
            grid on;
            hold off;
        end
        legend('Location', 'bestoutside');

        % 3. Relative Elements (Mean)
        figure;
        sgtitle(sprintf('Relative Orbital Elements (Mean) - %s', deputy_type));
        for idx = 1:length(rel_fields_mean)
            field_name = rel_fields_mean{idx};
            subplot(3,2,idx);
            hold on;
            curve_idx = 0;
            for init_idx = 1:length(cases_to_plot)
                for j2_idx = 1:length(j2_conditions)
                    curve_idx = curve_idx + 1;
                    init_case = cases_to_plot{init_idx};
                    j2_case = j2_conditions{j2_idx};
                    deputy_name = deputy_mapping.(init_case).(deputy_type);

                    values = roe_results.(init_case).(j2_case).(deputy_name).(field_name);

                    plot(time_orbits, values, 'DisplayName', custom_legends{curve_idx});
                end
            end
            %title(strrep(field_name, '_', '\_'));
            xlabel('Number of Orbits');
            ylabel(rel_labels_mean{idx});
            grid on;
            hold off;
        end
        legend('Location', 'bestoutside');

        % 4. Relative Elements (Oscillating)
        figure;
        sgtitle(sprintf('Relative Orbital Elements (Oscullating) - %s', deputy_type));
        for idx = 1:length(rel_fields_osc)
            field_name = rel_fields_osc{idx};
            subplot(3,2,idx);
            hold on;
            curve_idx = 0;
            for init_idx = 1:length(cases_to_plot)
                for j2_idx = 1:length(j2_conditions)
                    curve_idx = curve_idx + 1;
                    init_case = cases_to_plot{init_idx};
                    j2_case = j2_conditions{j2_idx};
                    deputy_name = deputy_mapping.(init_case).(deputy_type);

                    values = roe_results.(init_case).(j2_case).(deputy_name).(field_name);

                    plot(time_orbits, values, 'DisplayName', custom_legends{curve_idx});
                end
            end
            %title(strrep(field_name, '_', '\_'));
            xlabel('Number of Orbits');
            ylabel(rel_labels_osc{idx});
            grid on;
            hold off;
        end
        legend('Location', 'bestoutside');
    end

    % === RTN Projections for Each Deputy (using subplot, legend below) ===
    for deputy_idx = 1:length(all_deputies)
        deputy_type = all_deputies{deputy_idx};

        figure('Color', 'w');
        sgtitle(sprintf('RTN Relative Motion (Mean Elements) - %s', deputy_type));

        % Create subplots manually
        ax1 = subplot(1,3,1);
        hold(ax1, 'on');
        axis equal;
        grid on;
        xlabel('Tangential [m]');
        ylabel('Radial [m]');
        %title('R-T Plane');

        ax2 = subplot(1,3,2);
        hold(ax2, 'on');
        axis equal;
        grid on;
        xlabel('Normal [m]');
        ylabel('Radial [m]');
        %title('R-N Plane');

        ax3 = subplot(1,3,3);
        hold(ax3, 'on');
        axis equal;
        grid on;
        xlabel('Tangential [m]');
        ylabel('Normal [m]');
        %title('N-T Plane');

        % Preallocate handles for legend
        legend_handles = [];

        curve_idx = 0;
        for init_idx = 1:length(cases_to_plot)
            for j2_idx = 1:length(j2_conditions)
                curve_idx = curve_idx + 1;
                init_case = cases_to_plot{init_idx};
                j2_case = j2_conditions{j2_idx};
                deputy_name = deputy_mapping.(init_case).(deputy_type);

                % Extract rho_RTN
                rho_RTN = roe_results.(init_case).(j2_case).(deputy_name).rho_RTN;

                % Thicker lines for no J2
                if strcmp(j2_case, 'no_j2')
                    lw = 3; % thicker
                else
                    lw = 1.5; % thinner
                end

                % Plot into each subplot (convert km --> m)
                h = plot(ax1, rho_RTN(:,2)*1000, rho_RTN(:,1)*1000, 'LineWidth', lw); % T vs R
                plot(ax2, rho_RTN(:,3)*1000, rho_RTN(:,1)*1000, 'LineWidth', lw);       % N vs R
                plot(ax3, rho_RTN(:,2)*1000, rho_RTN(:,3)*1000, 'LineWidth', lw);       % T vs N

                % Save one handle per curve for legend
                legend_handles(end+1) = h;
            end
        end

        % Release holds
        hold(ax1, 'off');
        hold(ax2, 'off');
        hold(ax3, 'off');

        % Now manually create a legend underneath
        % Use normalized units to place the legend below all plots
        lgd = legend(legend_handles, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.3, 0.02, 0.4, 0.05]; % [left, bottom, width, height]
        lgd.FontSize = 10;
    end

    % === 3D RTN Relative Trajectories for Each Deputy ===
    for deputy_idx = 1:length(all_deputies)
        deputy_type = all_deputies{deputy_idx};

        figure('Color', 'w');
        hold on;
        grid on;
        axis equal;
        view(3); % 3D view
        xlabel('Radial [m]');
        ylabel('Tangential [m]');
        zlabel('Normal [m]');
        title(sprintf('3D Relative Motion (Mean Elements) - %s', deputy_type));
        
        % Preallocate handles for legend
        legend_handles = [];

        curve_idx = 0;
        for init_idx = 1:length(cases_to_plot)
            for j2_idx = 1:length(j2_conditions)
                curve_idx = curve_idx + 1;
                init_case = cases_to_plot{init_idx};
                j2_case = j2_conditions{j2_idx};
                deputy_name = deputy_mapping.(init_case).(deputy_type);

                % Extract rho_RTN
                rho_RTN = roe_results.(init_case).(j2_case).(deputy_name).rho_RTN;

                % Thicker lines for no-J2
                if strcmp(j2_case, 'no_j2')
                    lw = 2.5;
                else
                    lw = 1.5;
                end

                % 3D plot (km --> m)
                h = plot3(rho_RTN(:,1)*1e3, rho_RTN(:,2)*1e3, rho_RTN(:,3)*1e3, 'LineWidth', lw);
                legend_handles(end+1) = h;
            end
        end
        hold off;

        % Add legend
        lgd = legend(legend_handles, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.3, 0.02, 0.4, 0.05];
        lgd.FontSize = 10;
    end

        % === Relative Orbital Elements Projections for Each Deputy ===
    for deputy_idx = 1:length(all_deputies)
        deputy_type = all_deputies{deputy_idx};

        figure('Color', 'w');
        sgtitle(sprintf('Relative Orbital Elements (Osculating) - %s', deputy_type));

        % 1. Plot relative eccentricity vector (d_e_x vs d_e_y)
        ax1 = subplot(1,3,1);
        hold(ax1, 'on');
        axis equal;
        grid on;
        xlabel('\Delta e_x [m]');
        ylabel('\Delta e_y [m]');
        title('Relative Eccentricity Vector');

        % 2. Plot relative inclination vector (d_i_x vs d_i_y)
        ax2 = subplot(1,3,2);
        hold(ax2, 'on');
        axis equal;
        grid on;
        xlabel('\Delta i_x [m]');
        ylabel('\Delta i_y [m]');
        title('Relative Inclination Vector');

        % 3. Plot relative mean longitude vs semi-major axis (d_lambda vs d_a)
        ax3 = subplot(1,3,3);
        hold(ax3, 'on');
        axis equal;
        grid on;
        xlabel('\Delta \lambda [m]');
        ylabel('\Delta a [m]');
        title('Relative Mean Longitude vs Semi-Major Axis');

        % Preallocate handles for legend
        legend_handles_osc = [];

        curve_idx = 0;
        for init_idx = 1:length(cases_to_plot)
            for j2_idx = 1:length(j2_conditions)
                curve_idx = curve_idx + 1;
                init_case = cases_to_plot{init_idx};
                j2_case = j2_conditions{j2_idx};
                deputy_name = deputy_mapping.(init_case).(deputy_type);

                % Extract oscillating relative elements
                d_e_x_osc = roe_results.(init_case).(j2_case).(deputy_name).d_e_x_osc;
                d_e_y_osc = roe_results.(init_case).(j2_case).(deputy_name).d_e_y_osc;
                d_i_x_osc = roe_results.(init_case).(j2_case).(deputy_name).d_i_x_osc;
                d_i_y_osc = roe_results.(init_case).(j2_case).(deputy_name).d_i_y_osc;
                d_lambda_osc = roe_results.(init_case).(j2_case).(deputy_name).d_lambda_osc;
                d_a_osc = roe_results.(init_case).(j2_case).(deputy_name).d_a_osc;

                % Thicker lines for no-J2 cases
                if strcmp(j2_case, 'no_j2')
                    lw = 2.5;
                else
                    lw = 1.5;
                end

                % Plot each subplot
                h = plot(ax1, d_e_x_osc, d_e_y_osc, 'LineWidth', lw);
                plot(ax2, d_i_x_osc, d_i_y_osc, 'LineWidth', lw);
                plot(ax3, d_lambda_osc, d_a_osc, 'LineWidth', lw);

                % Save handle for legend
                legend_handles_osc(end+1) = h;
            end
        end

        hold(ax1, 'off');
        hold(ax2, 'off');
        hold(ax3, 'off');

        % Add legend below
        lgd = legend(legend_handles_osc, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.3, 0.02, 0.4, 0.05];
        lgd.FontSize = 10;

        % === Now repeat for Mean elements ===
        figure('Color', 'w');
        sgtitle(sprintf('Relative Orbital Elements (Mean) - %s', deputy_type));

        ax1 = subplot(1,3,1);
        hold(ax1, 'on');
        axis equal;
        grid on;
        xlabel('\Delta e_x [m]');
        ylabel('\Delta e_y [m]');
        title('Relative Eccentricity Vector (Mean)');

        ax2 = subplot(1,3,2);
        hold(ax2, 'on');
        axis equal;
        grid on;
        xlabel('\Delta i_x [m]');
        ylabel('\Delta i_y [m]');
        title('Relative Inclination Vector (Mean)');

        ax3 = subplot(1,3,3);
        hold(ax3, 'on');
        axis equal;
        grid on;
        xlabel('\Delta \lambda [m]');
        ylabel('\Delta a [m]');
        title('Relative Mean Longitude vs Semi-Major Axis (Mean)');

        % Preallocate handles
        legend_handles_mean = [];

        curve_idx = 0;
        for init_idx = 1:length(cases_to_plot)
            for j2_idx = 1:length(j2_conditions)
                curve_idx = curve_idx + 1;
                init_case = cases_to_plot{init_idx};
                j2_case = j2_conditions{j2_idx};
                deputy_name = deputy_mapping.(init_case).(deputy_type);

                % Extract mean relative elements
                d_e_x_mean = roe_results.(init_case).(j2_case).(deputy_name).d_e_x_mean;
                d_e_y_mean = roe_results.(init_case).(j2_case).(deputy_name).d_e_y_mean;
                d_i_x_mean = roe_results.(init_case).(j2_case).(deputy_name).d_i_x_mean;
                d_i_y_mean = roe_results.(init_case).(j2_case).(deputy_name).d_i_y_mean;
                d_lambda_mean = roe_results.(init_case).(j2_case).(deputy_name).d_lambda_mean;
                d_a_mean = roe_results.(init_case).(j2_case).(deputy_name).d_a_mean;

                if strcmp(j2_case, 'no_j2')
                    lw = 2.5;
                else
                    lw = 1.5;
                end

                % Plot
                h = plot(ax1, d_e_x_mean, d_e_y_mean, 'LineWidth', lw);
                plot(ax2, d_i_x_mean, d_i_y_mean, 'LineWidth', lw);
                plot(ax3, d_lambda_mean, d_a_mean, 'LineWidth', lw);

                % Save handle for legend
                legend_handles_mean(end+1) = h;
            end
        end

        hold(ax1, 'off');
        hold(ax2, 'off');
        hold(ax3, 'off');

        lgd = legend(legend_handles_mean, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.3, 0.02, 0.4, 0.05];
        lgd.FontSize = 10;
    end



end
