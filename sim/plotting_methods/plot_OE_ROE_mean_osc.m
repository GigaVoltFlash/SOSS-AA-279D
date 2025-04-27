function plot_OE_ROE_mean_osc(roe_results, tstart, tint, tend)
% plot_roe_results
% Inputs:
%   roe_results : struct containing ROE data (output from compute_OE_ROE_mean_osc and organized into struct)
%   tstart, tint, tend : simulation start time, timestep, and end time [s]

    % Build time vector
    time = tstart:tint:tend;
    
    % Cases and J2 conditions to plot
    cases_to_plot = {'initial1', 'initial2'};
    j2_conditions = {'no_j2', 'with_j2'};

    % Absolute orbital elements (quasi-nonsingular OE)
    abs_fields = {'a', 'e_x', 'e_y', 'i', 'RAAN', 'u'};

    % Relative orbital elements
    osc_fields = {'d_a_osc', 'd_lambda_osc', 'd_e_x_osc', 'd_e_y_osc', 'd_i_x_osc', 'd_i_y_osc'};
    mean_fields = {'d_a_mean', 'd_lambda_mean', 'd_e_x_mean', 'd_e_y_mean', 'd_i_x_mean', 'd_i_y_mean'};

    % === Subplots for Absolute Orbital Elements ===
    figure;
    sgtitle('Absolute Orbital Elements (Osculating and Mean)');
    for idx = 1:length(abs_fields)
        field_name = abs_fields{idx};
        subplot(3,2,idx);
        hold on;
        for init_idx = 1:length(cases_to_plot)
            for j2_idx = 1:length(j2_conditions)
                init_case = cases_to_plot{init_idx};
                j2_case = j2_conditions{j2_idx};

                % Mean and Oscillating
                mean_values = roe_results.(init_case).(j2_case).([field_name, '_mean']);
                osc_values  = roe_results.(init_case).(j2_case).([field_name, '_osc']);

                plot(time, mean_values, '--', 'DisplayName', sprintf('%s %s mean', init_case, j2_case));
                plot(time, osc_values, '-', 'DisplayName', sprintf('%s %s osc', init_case, j2_case));
            end
        end
        title(strrep(field_name, '_', '\_'));
        xlabel('Time [s]');
        ylabel(field_name);
        grid on;
        hold off;
    end
    legend('Location', 'bestoutside');
    
    % === Subplots for Relative Orbital Elements (Osculating) ===
    figure;
    sgtitle('Relative Orbital Elements (Osculating)');
    for idx = 1:length(osc_fields)
        field_name = osc_fields{idx};
        subplot(3,2,idx);
        hold on;
        for init_idx = 1:length(cases_to_plot)
            for j2_idx = 1:length(j2_conditions)
                init_case = cases_to_plot{init_idx};
                j2_case = j2_conditions{j2_idx};

                values = roe_results.(init_case).(j2_case).(field_name);

                plot(time, values, 'DisplayName', sprintf('%s %s', init_case, j2_case));
            end
        end
        title(strrep(field_name, '_', '\_'));
        xlabel('Time [s]');
        ylabel(field_name);
        grid on;
        hold off;
    end
    legend('Location', 'bestoutside');

    % === Subplots for Relative Orbital Elements (Mean) ===
    figure;
    sgtitle('Relative Orbital Elements (Mean)');
    for idx = 1:length(mean_fields)
        field_name = mean_fields{idx};
        subplot(3,2,idx);
        hold on;
        for init_idx = 1:length(cases_to_plot)
            for j2_idx = 1:length(j2_conditions)
                init_case = cases_to_plot{init_idx};
                j2_case = j2_conditions{j2_idx};

                values = roe_results.(init_case).(j2_case).(field_name);

                plot(time, values, 'DisplayName', sprintf('%s %s', init_case, j2_case));
            end
        end
        title(strrep(field_name, '_', '\_'));
        xlabel('Time [s]');
        ylabel(field_name);
        grid on;
        hold off;
    end
    legend('Location', 'bestoutside');
end
