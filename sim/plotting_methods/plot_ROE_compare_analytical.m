function plot_ROE_compare_analytical(roe_results, roe_analytical_results, tstart, tint, tend, t_orbit)
% Inputs:
%   roe_results : struct containing ROE data (organized with SV2 and SV3)
%   roe_analytical_results: struct containing ROE data from the analytical first order formulation.
%   tstart, tint, tend : simulation start time, timestep, and end time [s]

    % Build time vector
    time = tstart:tint:tend;
    time_orbits = time / t_orbit;

    % Legend entries you want
    % custom_legends = {'IC 1 J2 propagated', 'IC 1 J2 analytical', 'IC 2 J2 propagated','IC 2 J2 analytical'};
    custom_legends = {'IC 2 J2 propagated','IC 2 J2 analytical'};

    % Deputies
    all_deputies = {'SV2'};

    % Mapping of initial case to actual deputy names
    deputy_mapping = struct( ...
        'initial1', struct('SV2', 'SV2_1', 'SV3', 'SV3_1'), ...
        'initial2', struct('SV2', 'SV2_2', 'SV3', 'SV3_2'), ...
        'initial3', struct('SV2', 'SV2_3', 'SV3', 'SV3_3') ...
    );

    rel_fields = {'d_a_mean', 'd_lambda_mean', 'd_e_x_mean', 'd_e_y_mean', 'd_i_x_mean', 'd_i_y_mean'};
    rel_labels = {'\delta a [m]', '\delta \lambda [m]', '\delta e_x [m]', '\delta e_y [m]', '\delta i_x [m]', '\delta i_y [m]'};

    % === Now for each Deputy separately (SV2, SV3) ===
    for deputy_idx = 1:length(all_deputies)
        deputy_type = all_deputies{deputy_idx};

        % Running this once with the given ICs (over t)
        figure;
        sgtitle(sprintf('Relative Orbital Elements - %s', deputy_type));
        for idx = 1:length(rel_fields)
            field_name = rel_fields{idx};
            subplot(3,2,idx);
            hold on;
            deputy_name = deputy_mapping.initial2.(deputy_type);

            values = roe_results.initial2.with_j2.(deputy_name).(field_name);
            values_analytical_all = roe_analytical_results.initial2.(deputy_name).roe_analytical_j2_given_ic;
            values_analytical_field = values_analytical_all(:, idx);

            plot(time_orbits, values, 'DisplayName', custom_legends{1});
            plot(time_orbits, values_analytical_field, 'LineStyle', "--", 'DisplayName', custom_legends{2});
            xlabel('Number of Orbits');
            ylabel(rel_labels{idx});
            grid on;
            hold off;
        end
        legend();

        % saveas(gcf, sprintf('figures/PS4/ROE_time_analytical_compare_given_ICs_%s.png', deputy_type));

        figure;
        sgtitle(sprintf('Relative Orbital Elements - %s', deputy_type));
        for idx = 1:length(rel_fields)
            field_name = rel_fields{idx};
            subplot(3,2,idx);
            hold on;
            deputy_name = deputy_mapping.initial2.(deputy_type);

            values = roe_results.initial2.with_j2.(deputy_name).(field_name);
            values_analytical_all = roe_analytical_results.initial2.(deputy_name).roe_analytical_j2_mean_ic;
            values_analytical_field = values_analytical_all(:, idx);

            plot(time_orbits, values, 'DisplayName', custom_legends{1});
            plot(time_orbits, values_analytical_field, 'LineStyle', "--", 'DisplayName', custom_legends{2});
            xlabel('Number of Orbits');
            ylabel(rel_labels{idx});
            grid on;
            hold off;
        end
        legend();

        % saveas(gcf, sprintf('figures/PS4/ROE_time_analytical_compare_mean_ICs_%s.png', deputy_type));

        % Plot the ROE space plots that compare the propagated and the 
        plot_ROE_space_compare_analytical(roe_results, roe_analytical_results, deputy_mapping, deputy_type, true, 2);
        plot_ROE_space_compare_analytical(roe_results, roe_analytical_results, deputy_mapping, deputy_type, false, 2);
        plot_ROE_space_compare_analytical(roe_results, roe_analytical_results, deputy_mapping, deputy_type, true, 3);
        plot_ROE_space_compare_analytical(roe_results, roe_analytical_results, deputy_mapping, deputy_type, false, 3);
    end
end