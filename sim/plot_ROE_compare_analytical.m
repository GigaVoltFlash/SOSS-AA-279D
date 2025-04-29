function plot_ROE_compare_analytical(roe_results, roe_analytical_results, tstart, tint, tend, t_orbit)
% plot_OE_ROE_mean_osc
% Inputs:
%   roe_results : struct containing ROE data (organized with SV2 and SV3)
%   tstart, tint, tend : simulation start time, timestep, and end time [s]

    % Build time vector
    time = tstart:tint:tend;
    time_orbits = time / t_orbit;
    
    % Cases and J2 conditions
    % cases_to_plot = {'initial1', 'initial2'};
    cases_to_plot = {'initial2'};

    % Legend entries you want
    % custom_legends = {'IC 1 J2 propagated', 'IC 1 J2 analytical', 'IC 2 J2 propagated','IC 2 J2 analytical'};
    custom_legends = {'IC 2 J2 propagated','IC 2 J2 analytical'};

    % Deputies
    all_deputies = {'SV2', 'SV3'};

    % Mapping of initial case to actual deputy names
    deputy_mapping = struct( ...
        'initial1', struct('SV2', 'SV2_1', 'SV3', 'SV3_1'), ...
        'initial2', struct('SV2', 'SV2_2', 'SV3', 'SV3_2') ...
    );

    rel_fields = {'d_a_mean', 'd_lambda_mean', 'd_e_x_mean', 'd_e_y_mean', 'd_i_x_mean', 'd_i_y_mean'};
    rel_labels = {'\Delta a [m]', '\Delta \lambda [m]', '\Delta e_x [m]', '\Delta e_y [m]', '\Delta i_x [m]', '\Delta i_y [m]'};

    % === Now for each Deputy separately (SV2, SV3) ===
    for deputy_idx = 1:length(all_deputies)
        deputy_type = all_deputies{deputy_idx};

        figure;
        sgtitle(sprintf('Relative Orbital Elements - %s', deputy_type));
        for idx = 1:length(rel_fields)
            field_name = rel_fields{idx};
            subplot(3,2,idx);
            hold on;
            curve_idx = 1;
            for init_idx = 1:length(cases_to_plot)
                init_case = cases_to_plot{init_idx};
                deputy_name = deputy_mapping.(init_case).(deputy_type);

                values = roe_results.(init_case).with_j2.(deputy_name).(field_name);
                values_analytical_all = roe_analytical_results.(init_case).(deputy_name).roe_analytical_j2;
                values_analytical_field = values_analytical_all(:, idx);

                plot(time_orbits, values, 'DisplayName', custom_legends{curve_idx});
                plot(time_orbits, values_analytical_field, 'LineStyle', "--", 'DisplayName', custom_legends{curve_idx+1});
                curve_idx = curve_idx + 2;
            end
            xlabel('Number of Orbits');
            ylabel(rel_labels{idx});
            grid on;
            hold off;
        end
        legend('Location', 'bestoutside');
    end
end
