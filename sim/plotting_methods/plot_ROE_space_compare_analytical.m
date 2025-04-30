% Helper function specifically made for use in plot_ROE_compare_analytical.
% Refer to that file.
function plot_ROE_space_compare_analytical(roe_results, roe_analytical_results, deputy_mapping, deputy_type, is_mean, ic_number)

    figure('Color', 'w','Position', [100, 100, 1000, 600]);
    ax1 = subplot(1,3,1);
    hold(ax1, 'on');
    axis equal;
    grid on;
    xlabel('\delta e_x [m]');
    ylabel('\delta e_y [m]');
    % title('Relative Eccentricity Vector (Mean)');

    ax2 = subplot(1,3,2);
    hold(ax2, 'on');
    axis equal;
    grid on;
    xlabel('\delta i_x [m]');
    ylabel('\delta i_y [m]');
    % title('Relative Inclination Vector (Mean)');

    ax3 = subplot(1,3,3);
    hold(ax3, 'on');
    axis equal;
    grid on;
    xlabel('\delta \lambda [m]');
    ylabel('\delta a [m]');
    % title('Relative Mean Longitude vs Semi-Major Axis (Mean)');

    % Extract mean relative elements
    if ic_number == 2
        ic = 'initial2';
    else
        ic = 'initial3';
    end
    deputy_name = deputy_mapping.(ic).(deputy_type);


    d_e_x_mean = roe_results.(ic).with_j2.(deputy_name).d_e_x_mean;
    d_e_y_mean = roe_results.(ic).with_j2.(deputy_name).d_e_y_mean;
    d_i_x_mean = roe_results.(ic).with_j2.(deputy_name).d_i_x_mean;
    d_i_y_mean = roe_results.(ic).with_j2.(deputy_name).d_i_y_mean;
    d_lambda_mean = roe_results.(ic).with_j2.(deputy_name).d_lambda_mean;
    d_a_mean = roe_results.(ic).with_j2.(deputy_name).d_a_mean;

    % Extract mean relative elements
    if is_mean
        values_analytical_all = roe_analytical_results.(ic).(deputy_name).roe_analytical_j2_mean_ic;
    else
        values_analytical_all = roe_analytical_results.(ic).(deputy_name).roe_analytical_j2_given_ic;
    end
    d_e_x_analytical = values_analytical_all(:,3);
    d_e_y_analytical = values_analytical_all(:,4);
    d_i_x_analytical = values_analytical_all(:,5);
    d_i_y_analytical = values_analytical_all(:,6);
    d_lambda_analytical = values_analytical_all(:,2);
    d_a_analytical = values_analytical_all(:,1);

    % Plot
    h = plot(ax1, d_e_x_mean, d_e_y_mean, 'DisplayName', 'Propagated', 'LineWidth', 1.5);
    plot(ax1, d_e_x_analytical, d_e_y_analytical, 'LineStyle', '--', 'DisplayName', 'Analytical', 'LineWidth', 1.5);
    plot(ax2, d_i_x_mean, d_i_y_mean, 'DisplayName', 'Propagated', 'LineWidth', 1.5);
    plot(ax2, d_i_x_analytical, d_i_y_analytical,'LineStyle', '--', 'DisplayName', 'Analytical', 'LineWidth', 1.5);
    plot(ax3, d_lambda_mean, d_a_mean, 'DisplayName', 'Propagated', 'LineWidth', 1.5);
    plot(ax3, d_lambda_analytical, d_a_analytical, 'LineStyle', '--', 'DisplayName', 'Analytical', 'LineWidth', 1.5);

    hold(ax1, 'off');
    hold(ax2, 'off');
    hold(ax3, 'off');

    lgd = legend('Orientation', 'horizontal');
    lgd.Units = 'normalized';
    lgd.Position = [0.3, 0.0, 0.4, 0.05];
    lgd.FontSize = 10;

    if is_mean
        saveas(gcf, sprintf('figures/PS4/ROE_analytical_compare_mean_IC%d_%s.png', ic_number, deputy_type));
    else
        saveas(gcf, sprintf('figures/PS4/ROE_analytical_compare_given_IC%d_%s.png', ic_number, deputy_type));
    end
end