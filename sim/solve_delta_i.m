function delta_i = solve_delta_i(a, delta_r_min, delta_e)
    % Solves for the magnitude of delta i that gives a desired minimum relative distance
    % Inputs:
    %   a            = semi-major axis (meters)
    %   delta_r_min = desired minimum distance (meters)
    %   delta_e     = chosen scalar value for |delta e|
    % Output:
    %   delta_i     = solution for |delta i| (if found)

    % Define the equation based on Eq. (2.22)
    f = @(di) (sqrt(2) * a * abs(delta_e * di)) / ...
              sqrt(delta_e^2 + di^2 + abs(delta_e + di) * abs(delta_e - di)) - delta_r_min;

    % Use a reasonable initial guess and bounds
    try
        delta_i = fzero(f, [1e-6, 0.01]);  % bounds from µrad to mrad
    catch
        warning('No solution found for given parameters.');
        delta_i = NaN;
    end
end
