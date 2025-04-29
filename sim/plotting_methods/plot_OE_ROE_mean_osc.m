function plot_OE_ROE_mean_osc(roe_results, tstart, tint, tend, t_orbit)
% Inputs:
%   roe_results : struct containing ROE data (organized with SV2 and SV3)
%   tstart, tint, tend : simulation start time, timestep, and end time [s]

    time = tstart:tint:tend;
    time_orbits = time / t_orbit;
    
    cases_to_plot = {'initial2'}; % , 'initial1'
    j2_conditions = {'no_j2', 'with_j2'};

    custom_legends = {'without J2', 'with J2'}; % 'IC 1 without J2', 'IC 1 with J2',

    all_deputies = {'SV2'};%'SV3'

    deputy_mapping = struct( ...
        'initial1', struct('SV2', 'SV2_1', 'SV3', 'SV3_1'), ...
        'initial2', struct('SV2', 'SV2_2', 'SV3', 'SV3_2') ...
    );

    abs_fields = {'a', 'e_x', 'e_y', 'i', 'RAAN', 'u'};
    abs_labels = {'a [km]', 'e_x', 'e_y', 'i [deg]', 'RAAN [deg]', 'u [deg]'};

    rel_fields_osc = {'d_a_osc', 'd_lambda_osc', 'd_e_x_osc', 'd_e_y_osc', 'd_i_x_osc', 'd_i_y_osc'};
    rel_labels_osc = {'\Delta a [m]', '\Delta \lambda [m]', '\Delta e_x [m]', '\Delta e_y [m]', '\Delta i_x [m]', '\Delta i_y [m]'};

    rel_fields_mean = {'d_a_mean', 'd_lambda_mean', 'd_e_x_mean', 'd_e_y_mean', 'd_i_x_mean', 'd_i_y_mean'};
    rel_labels_mean = {'\Delta a [m]', '\Delta \lambda [m]', '\Delta e_x [m]', '\Delta e_y [m]', '\Delta i_x [m]', '\Delta i_y [m]'};

    % Plot for each deputy
    for deputy_idx = 1:length(all_deputies)
        deputy_type = all_deputies{deputy_idx};

        % 1. Absolute Elements (Mean)
        figure('Position', [100, 100, 800, 600]);
        %sgtitle(sprintf('Absolute Orbital Elements (Mean) - %s', deputy_type));
        
        legend_handles = [];

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

                    values = roe_results.(init_case).(j2_case).(deputy_name).([field_name, '_mean']);
                    
                    if idx == 1
                        h = plot(time_orbits, values, 'LineWidth', 1.5);
                        legend_handles(end+1) = h;
                    else
                        plot(time_orbits, values, 'LineWidth', 1.5);
                    end
                end
            end
            xlabel('Number of Orbits');
            ylabel(abs_labels{idx});
            ylim_current = ylim;
            y_range = ylim_current(2) - ylim_current(1); 
            margin = 0.05 * y_range; 
            ylim([ylim_current(1) - margin, ylim_current(2) + margin]); 
            grid on;
            hold off;
        end
        lgd = legend(legend_handles, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.35, 0.02, 0.3, 0.03]; 
        lgd.FontSize = 10;

        saveas(gcf, sprintf('figures/PS4/OE_abs_mean_%s.png', deputy_type));

        % 2. Absolute Elements (Oscillating)
        figure('Position', [100, 100, 800, 600]);
        %sgtitle(sprintf('Absolute Orbital Elements (Oscullating) - %s', deputy_type));
        
        legend_handles = [];
        
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

                    if idx == 1
                        h = plot(time_orbits, values, 'LineWidth', 1.5);
                        legend_handles(end+1) = h;
                    else
                        plot(time_orbits, values, 'LineWidth', 1.5);
                    end
                end
            end
            xlabel('Number of Orbits');
            ylabel(abs_labels{idx});
            ylim_current = ylim;
            y_range = ylim_current(2) - ylim_current(1); 
            margin = 0.05 * y_range; 
            ylim([ylim_current(1) - margin, ylim_current(2) + margin]); 
            grid on;
            hold off;
        end
        lgd = legend(legend_handles, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.35, 0.02, 0.3, 0.03]; 
        lgd.FontSize = 10;

        saveas(gcf, sprintf('figures/PS4/OE_abs_osc_%s.png', deputy_type));


        % 3. Relative Elements (Mean)
        figure('Position', [100, 100, 800, 600]);
        %sgtitle(sprintf('Relative Orbital Elements (Mean) - %s', deputy_type));

        legend_handles = [];

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

                    if idx == 1
                        h = plot(time_orbits, values, 'LineWidth', 1.5);
                        legend_handles(end+1) = h;
                    else
                        plot(time_orbits, values, 'LineWidth', 1.5);
                    end
                end
            end
            xlabel('Number of Orbits');
            ylabel(rel_labels_mean{idx});
            ylim_current = ylim;
            y_range = ylim_current(2) - ylim_current(1); 
            margin = 0.05 * y_range; 
            ylim([ylim_current(1) - margin, ylim_current(2) + margin]); 
            grid on;
            hold off;
        end
        lgd = legend(legend_handles, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.35, 0.02, 0.3, 0.03]; 
        lgd.FontSize = 10;

        saveas(gcf, sprintf('figures/PS4/ROE_mean_%s.png', deputy_type));

        % 4. Relative Elements (Oscillating)
        figure('Position', [100, 100, 800, 600]);
        %sgtitle(sprintf('Relative Orbital Elements (Oscullating) - %s', deputy_type));

        legend_handles = [];

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

                    if idx == 1
                        h = plot(time_orbits, values, 'LineWidth', 1.5);
                        legend_handles(end+1) = h;
                    else
                        plot(time_orbits, values, 'LineWidth', 1.5);
                    end
                end
            end
            xlabel('Number of Orbits');
            ylabel(rel_labels_osc{idx});
            ylim_current = ylim;
            y_range = ylim_current(2) - ylim_current(1); 
            margin = 0.05 * y_range; 
            ylim([ylim_current(1) - margin, ylim_current(2) + margin]); 
            grid on;
            hold off;
        end
        lgd = legend(legend_handles, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.35, 0.02, 0.3, 0.03]; 
        lgd.FontSize = 10;

        saveas(gcf, sprintf('figures/PS4/ROE_osc_%s.png', deputy_type));
    end

    % RTN Projections using mean
    for deputy_idx = 1:length(all_deputies)
        deputy_type = all_deputies{deputy_idx};

        figure('Color', 'w','Position', [100, 100, 800, 600]);
        %sgtitle(sprintf('RTN Relative Motion (Mean Elements) - %s', deputy_type));

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

        legend_handles = [];

        curve_idx = 0;
        for init_idx = 1:length(cases_to_plot)
            for j2_idx = 1:length(j2_conditions)
                curve_idx = curve_idx + 1;
                init_case = cases_to_plot{init_idx};
                j2_case = j2_conditions{j2_idx};
                deputy_name = deputy_mapping.(init_case).(deputy_type);

                rho_RTN = roe_results.(init_case).(j2_case).(deputy_name).rho_RTN;

                if strcmp(j2_case, 'no_j2')
                    lw = 2.5;
                    line_style = '-';
                else
                    lw = 1.5;
                    line_style = '-';
                end

                % Plot into each subplot (convert km --> m)
                h = plot(ax1, rho_RTN(:,2)*1000, rho_RTN(:,1)*1000, line_style,'LineWidth', lw); % T vs R
                plot(ax2, rho_RTN(:,3)*1000, rho_RTN(:,1)*1000, line_style,'LineWidth', lw);       % N vs R
                plot(ax3, rho_RTN(:,2)*1000, rho_RTN(:,3)*1000, line_style,'LineWidth', lw);       % T vs N

                legend_handles(end+1) = h;
            end
        end

        hold(ax1, 'off');
        hold(ax2, 'off');
        hold(ax3, 'off');

        lgd = legend(legend_handles, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.35, 0.02, 0.3, 0.03]; 
        lgd.FontSize = 10;

        saveas(gcf, sprintf('figures/PS4/RTN_projections_%s.png', deputy_type));
    end

    % 3D RTN
    for deputy_idx = 1:length(all_deputies)
        deputy_type = all_deputies{deputy_idx};

        figure('Color', 'w','Position', [100, 100, 800, 600]);
        hold on;
        grid on;
        axis equal;
        view(3); 
        xlabel('Radial [m]');
        ylabel('Tangential [m]');
        zlabel('Normal [m]');
        %title(sprintf('3D Relative Motion (Mean Elements) - %s', deputy_type));
        
        legend_handles = [];

        curve_idx = 0;
        for init_idx = 1:length(cases_to_plot)
            for j2_idx = 1:length(j2_conditions)
                curve_idx = curve_idx + 1;
                init_case = cases_to_plot{init_idx};
                j2_case = j2_conditions{j2_idx};
                deputy_name = deputy_mapping.(init_case).(deputy_type);

                rho_RTN = roe_results.(init_case).(j2_case).(deputy_name).rho_RTN;

                if strcmp(j2_case, 'no_j2')
                    lw = 2.5;
                    line_style='-';
                else
                    lw = 1.5;
                    line_style='-';
                end

                % 3D plot (km --> m)
                h = plot3(rho_RTN(:,1)*1e3, rho_RTN(:,2)*1e3, rho_RTN(:,3)*1e3,line_style, 'LineWidth', lw);
                legend_handles(end+1) = h;
            end
        end
        hold off;

        lgd = legend(legend_handles, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.35, 0.02, 0.3, 0.03]; 
        lgd.FontSize = 10;

        saveas(gcf, sprintf('figures/PS4/RTN_3D_%s.png', deputy_type));
    end

    % ROE Projections
    for deputy_idx = 1:length(all_deputies)
        deputy_type = all_deputies{deputy_idx};

        % Osculating
        figure('Color', 'w','Position', [100, 100, 1000, 600]);
        %sgtitle(sprintf('Relative Orbital Elements (Osculating) - %s', deputy_type));

        % Plot relative eccentricity vector (d_e_x vs d_e_y)
        ax1 = subplot(1,3,1);
        hold(ax1, 'on');
        axis equal;
        grid on;
        xlabel('\Delta e_x [m]');
        ylabel('\Delta e_y [m]');
        %title('Relative Eccentricity Vector');

        % Plot relative inclination vector (d_i_x vs d_i_y)
        ax2 = subplot(1,3,2);
        hold(ax2, 'on');
        axis equal;
        grid on;
        xlabel('\Delta i_x [m]');
        ylabel('\Delta i_y [m]');
        %title('Relative Inclination Vector');

        % Plot relative mean longitude vs semi-major axis (d_lambda vs d_a)
        ax3 = subplot(1,3,3);
        hold(ax3, 'on');
        axis equal;
        grid on;
        xlabel('\Delta \lambda [m]');
        ylabel('\Delta a [m]');
        %title('Relative Mean Longitude vs Semi-Major Axis');

        legend_handles_osc = [];

        curve_idx = 0;
        for init_idx = 1:length(cases_to_plot)
            for j2_idx = 1:length(j2_conditions)
                curve_idx = curve_idx + 1;
                init_case = cases_to_plot{init_idx};
                j2_case = j2_conditions{j2_idx};
                deputy_name = deputy_mapping.(init_case).(deputy_type);

                d_e_x_osc = roe_results.(init_case).(j2_case).(deputy_name).d_e_x_osc;
                d_e_y_osc = roe_results.(init_case).(j2_case).(deputy_name).d_e_y_osc;
                d_i_x_osc = roe_results.(init_case).(j2_case).(deputy_name).d_i_x_osc;
                d_i_y_osc = roe_results.(init_case).(j2_case).(deputy_name).d_i_y_osc;
                d_lambda_osc = roe_results.(init_case).(j2_case).(deputy_name).d_lambda_osc;
                d_a_osc = roe_results.(init_case).(j2_case).(deputy_name).d_a_osc;

                if strcmp(j2_case, 'no_j2')
                    lw = 2.5;
                    
                    h = plot(ax1, d_e_x_osc, d_e_y_osc, 'LineWidth', lw);
                    plot(ax2, d_i_x_osc, d_i_y_osc, 'LineWidth', lw);
                    plot(ax3, d_lambda_osc, d_a_osc, 'LineWidth', lw);

                    plot(ax1, d_e_x_osc(1), d_e_y_osc(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', h.Color, 'MarkerEdgeColor', h.Color);
                    plot(ax2, d_i_x_osc(1), d_i_y_osc(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', h.Color, 'MarkerEdgeColor', h.Color);
                    plot(ax3, d_lambda_osc(1), d_a_osc(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', h.Color, 'MarkerEdgeColor', h.Color);
                else
                    lw = 1.5;

                    h = plot(ax1, d_e_x_osc, d_e_y_osc,'Color',"#D95319", 'LineWidth', lw);
                    plot(ax2, d_i_x_osc, d_i_y_osc,'Color', "#D95319", 'LineWidth', lw);
                    plot(ax3, d_lambda_osc, d_a_osc,'Color', "#D95319", 'LineWidth', lw);
                end

                legend_handles_osc(end+1) = h;
            end
        end

        hold(ax1, 'off');
        hold(ax2, 'off');
        hold(ax3, 'off');

        lgd = legend(legend_handles_osc, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.35, 0.02, 0.3, 0.03]; 
        lgd.FontSize = 10;

        saveas(gcf, sprintf('figures/PS4/ROE_projections_osc_%s.png', deputy_type));

        % Mean
        figure('Color', 'w','Position', [100, 100, 1000, 600]);
        %sgtitle(sprintf('Relative Orbital Elements (Mean) - %s', deputy_type));

        ax1 = subplot(1,3,1);
        hold(ax1, 'on');
        axis equal;
        grid on;
        xlabel('\Delta e_x [m]');
        ylabel('\Delta e_y [m]');
        %title('Relative Eccentricity Vector (Mean)');

        ax2 = subplot(1,3,2);
        hold(ax2, 'on');
        axis equal;
        grid on;
        xlabel('\Delta i_x [m]');
        ylabel('\Delta i_y [m]');
        %title('Relative Inclination Vector (Mean)');

        ax3 = subplot(1,3,3);
        hold(ax3, 'on');
        axis equal;
        grid on;
        xlabel('\Delta \lambda [m]');
        ylabel('\Delta a [m]');
        %title('Relative Mean Longitude vs Semi-Major Axis (Mean)');

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
                    
                    h = plot(ax1, d_e_x_mean, d_e_y_mean, 'LineWidth', lw);
                    plot(ax2, d_i_x_mean, d_i_y_mean, 'LineWidth', lw);
                    plot(ax3, d_lambda_mean, d_a_mean, 'LineWidth', lw);

                    plot(ax1, d_e_x_mean(1), d_e_y_mean(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', h.Color, 'MarkerEdgeColor', h.Color);
                    plot(ax2, d_i_x_mean(1), d_i_y_mean(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', h.Color, 'MarkerEdgeColor', h.Color);
                    plot(ax3, d_lambda_mean(1), d_a_mean(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', h.Color, 'MarkerEdgeColor', h.Color);
                else
                    lw = 1.5;

                    h = plot(ax1, d_e_x_mean, d_e_y_mean,'Color',"#D95319", 'LineWidth', lw);
                    plot(ax2, d_i_x_mean, d_i_y_mean,'Color', "#D95319", 'LineWidth', lw);
                    plot(ax3, d_lambda_mean, d_a_mean,'Color', "#D95319", 'LineWidth', lw);
                end

                % Save handle for legend
                legend_handles_mean(end+1) = h;
            end
        end

        hold(ax1, 'off');
        hold(ax2, 'off');
        hold(ax3, 'off');

        lgd = legend(legend_handles_mean, custom_legends, 'Orientation', 'horizontal');
        lgd.Units = 'normalized';
        lgd.Position = [0.35, 0.02, 0.3, 0.03]; 
        lgd.FontSize = 10;

        saveas(gcf, sprintf('figures/PS4/ROE_projections_mean_%s.png', deputy_type));
    end
end
