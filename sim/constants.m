%% AA279D Constants
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

function constants()
    global R_earth DU_earth AU mu_earth mu_mars mu_sun R_sun
    R_earth = 6378.13; % km
    DU_earth = R_earth; % km
    AU = 1.496e8; % km
    mu_earth = 3.986e5; % km^3/s^2
    mu_mars = 4.283e4; % km^3/s^2
    mu_sun = 1.327e11; % km^3/s^2
    R_sun = 6.963e5; % km
end
