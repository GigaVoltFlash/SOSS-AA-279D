function rho = atmospheric_density_Vallado(alt_km)
% ATMOSPHERIC_DENSITY_VALLADO Interpolates atmospheric density from Vallado Table 8-4
%
% Inputs:
%   alt_km - Altitude in km
%
% Outputs:
%   rho - Atmospheric density in kg/m^3

    % Vallado Table 8-4 data
    alt_table = [...
        150, 200, 250, 300, 350, 400, 450, 500, ...
        600, 700, 800, 900, 1000];

    rho_table = [...
        1.96e-9, 5.57e-10, 2.53e-10, 1.17e-10, 5.54e-11, 2.79e-11, ...
        1.50e-11, 7.52e-12, 2.22e-12, 9.37e-13, 4.10e-13, 1.89e-13, ...
        9.51e-14];

    % Extrapolate below 150 km and above 1000 km
    if alt_km <= alt_table(1)
        rho = rho_table(1);
    elseif alt_km >= alt_table(end)
        rho = rho_table(end);
    else
        % Linear interpolation
        rho = interp1(alt_table, rho_table, alt_km, 'linear');
    end
    
end