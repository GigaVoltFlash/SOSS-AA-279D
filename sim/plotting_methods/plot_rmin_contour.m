function plot_rmin_contour(a, level_set_values)
    % Inputs:
    %   a = semi-major axis in meters
    %   level_set_values = vector of r_min level set values to isolate in 2D

    % Define ranges
    delta_e_vals = linspace(1e-7, 2e-4, 300);  % 0.00001 to 0.001
    delta_i_vals = linspace(1e-7, 2e-4, 300);

    % Meshgrid
    [DE, DI] = meshgrid(delta_e_vals, delta_i_vals);

    % Compute r_min surface
    Numerator = sqrt(2) * a .* abs(DE .* DI);
    Denominator = sqrt(DE.^2 + DI.^2 + abs(DE + DI) .* abs(DE - DI));
    R_min = Numerator ./ Denominator;

    % 3D surface plot
    figure;
    surf(DE, DI, R_min, 'EdgeColor', 'none');
    xlabel('|\delta e|');
    ylabel('|\delta i|');
    zlabel('\delta r_{min} [m]');
    title(sprintf('\\delta r_{min} Surface for a = %.0f km', a/1e3));
    colorbar;
    grid on;
    view(45, 30);
    hold on;

    % Projected contours
    contour3(DE, DI, R_min, level_set_values, 'k', 'LineWidth', 1.2);

    % Transparent level set planes
    cmap = lines(length(level_set_values));
    for i = 1:length(level_set_values)
        z_level = level_set_values(i);
        fill3([min(delta_e_vals), max(delta_e_vals), max(delta_e_vals), min(delta_e_vals)], ...
              [min(delta_i_vals), min(delta_i_vals), max(delta_i_vals), max(delta_i_vals)], ...
              z_level * ones(1,4), cmap(i,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    hold off;
    legend('r_{min} Surface', 'Level Set Contours', 'Location', 'northeast');
    saveas(gcf, 'figures/PS5/min_dist_contour.png');

    % Separate 2D level set plots
    % for i = 1:length(level_set_values)
    %     figure;
    %     contour(a*DE, a*DI, R_min, [level_set_values(i) level_set_values(i)], 'LineWidth', 2);
    %     xlabel('a|\delta e|');
    %     ylabel('a|\delta i|');
    %     title(sprintf('Level Set: \\delta r_{min} = %.1f m at a = %.0f km', ...
    %           level_set_values(i), a/1e3));
    %     grid on;
    %     axis equal;
    % end

    
end
