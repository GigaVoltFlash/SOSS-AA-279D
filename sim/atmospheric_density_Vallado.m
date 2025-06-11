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
        2.070e-9, 2.789e-10, 7.248e-11, 2.418e-11, 9.518e-12, 3.725e-12, ...
        1.585e-12, 6.967e-13, 1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15, ...
        3.019e-15];

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