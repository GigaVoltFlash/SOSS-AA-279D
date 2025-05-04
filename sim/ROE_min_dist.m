function design_safe_ROE(delta_r_min, a)
    % DESIGN_SAFE_ROE finds valid (de, di) pairs that meet a passive safety minimum separation
    % Inputs:
    %   delta_r_min   = desired minimum separation in meters
    %   altitude_km   = chief orbit altitude in km (assumes circular orbit around Earth)
    
    % Constants
    %Re = 6378e3;             % Earth radius in meters
    %a = Re + altitude_km*1e3; % Semi-major axis in meters

    % Range of delta e values to test
    de_vals = linspace(1e-5, 0.005, 1000); % try 10 µrad to 10 mrad
    di_vals = NaN(size(de_vals));         % to store matching delta i values

    % Search for di such that passive safety constraint is met
    for i = 1:length(de_vals)
        de = de_vals(i);
        
        % Solve Eq. (2.22) for di numerically given de and desired delta_r_min
        f = @(di) (sqrt(2)*a*abs(de*di) / sqrt(de^2 + di^2 + abs(de + di)*abs(de - di))) - delta_r_min;
        try
            % Use initial bounds for root finding
            di = fzero(f, [1e-5, 0.01]);
            di_vals(i) = di;
        catch
            % Skip if no root found
        end
    end

    % Clean up
    valid = ~isnan(di_vals);
    de_plot = de_vals(valid);
    di_plot = di_vals(valid);
    
    % Plot results
    figure;
    plot(de_plot, di_plot, 'b', 'LineWidth', 2);
    xlabel('|\delta e|');
    ylabel('|\delta i|');
    %title(sprintf('ROE pairs for \\delta r_{min} ≥ %.1f m at %.0f km altitude', delta_r_min, altitude_km));
    grid on;
    axis equal;
end
