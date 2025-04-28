function plot_OE_ROE_mean_osc(roe_results, tstart, tint, tend, t_orbit)
% plot_OE_ROE_mean_osc
% Inputs:
%   roe_results : struct containing ROE data (organized with SV2 and SV3)
%   tstart, tint, tend : simulation start time, timestep, and end time [s]

    % Build time vector
    time = tstart:tint:tend;
    time_orbits = time / t_orbit;
    
    % Cases and J2 conditions
    cases_to_plot = {'initial1', 'initial2'};
    j2_conditions = {'no_j2', 'with_j2'};

    % Legend entries you want
    custom_legends = {'IC 1 without J2', 'IC 1 with J2', 'IC 2 without J2', 'IC 2 with J2'};

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
end
