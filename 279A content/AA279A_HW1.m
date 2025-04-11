%% Space Mechanics Homework 1
% Tycho Bogdanowitsch
clc; clear;
%% Constants
G = 6.67e-20; % km^3/kg/s^2
AU = 149597870.7; % km
m_sun = 1.99e30; % kg
m_moon = 7.35e22; % kg
m_planets = [0.33e24,4.87e24,0.642e24,...
    1898e24,568e24,86.8e24,102e24]; % kg
m_earth = 5.97e24; % kg
r_sun_earth = 1*AU; % km
r_moon_earth = 3.84e5; % km
r_sun_planets = AU*[0.39, 0.72, 1.52, 5.2, 9.54, 19.19, 30.06]; % km

mu_earth = 3.986e5; % km^3/s^2
mu_sun = 1.327e11; % km^3/s^2
R_earth = 6371; % km

%% Problem 2
disp('')
disp('Problem 2')
% Part 1
altitudes = [600,12000,35786]; % km

r_satellite_earth = R_earth + altitudes; % km

calc_accel = @(m_body, r_body) ...
    abs(G * m_body * ((1 ./(r_body - r_sun_earth - r_satellite_earth).^2) - ...
    (1 ./(r_body - r_sun_earth).^2)));

acceleration_earth = G * m_earth ./ (r_satellite_earth).^2;
acceleration_sun = G * m_sun * ((1 ./(r_sun_earth - r_satellite_earth).^2) - ...
    (1 ./(r_sun_earth).^2));
acceleration_moon = G * m_moon * ((1 ./(r_moon_earth - r_satellite_earth).^2) - ...
    (1 ./(r_moon_earth).^2));
acceleration_planets = zeros(length(m_planets),length(altitudes));

for i = 1:length(m_planets)
    acceleration_planets(i,:) = calc_accel(m_planets(i),r_sun_planets(i));
end

perturbations = [acceleration_earth;acceleration_sun;acceleration_moon;acceleration_planets];
names = [{'Earth'},{'Sun'},{'Moon'},{'Mercury'},{'Venus'},{'Mars'},...
    {'Jupiter'},{'Saturn'},{'Uranus'},{'Neptune'}];
perturbations_table = array2table(perturbations,...
    'VariableNames',{'LEO (600 km)','MEO (12,000 km)','GEO (35,786 km)'},...
    'RowNames',names);

disp('Perturbation Accelerations (km/s^2):');
disp(perturbations_table)

% Part 2
total_accel_LEO = acceleration_sun(1) + acceleration_moon(1)...
    + sum(acceleration_planets(:,1));
total_accel_GEO = acceleration_sun(3) + acceleration_moon(3)...
    + sum(acceleration_planets(:,3));

fprintf(['Total perturbation acceleration at 600 km altitude (LEO): ' ...
    '%.3e km/s^2\n'], total_accel_LEO);
fprintf(['Total perturbation acceleration at 35,786 km altitude (GEO): ' ...
    '%.3e km/s^2\n'], total_accel_GEO);
fprintf(['The total perturbation acceleration at GEO is %.3f times larger ' ...
    'than at LEO\n'], total_accel_GEO/total_accel_LEO);

%% Problem 3
disp('')
disp('Problem 3')
% Part 1
h_NS_apogee = 100; % km
e_NS = 0.87; 

r_NS_apogee = R_earth + h_NS_apogee; 
a_NS = r_NS_apogee / (1 + e_NS); 

v_NS_apogee = sqrt(mu_earth * (2 / r_NS_apogee - 1 / a_NS));

v_NS_circular = sqrt(mu_earth / r_NS_apogee);

fprintf('New Shepard Velocity at Apogee: %.3f km/s\n', v_NS_apogee);
fprintf('Required Circular Orbit Velocity: %.3f km/s\n', v_NS_circular);
fprintf(['The New Shepard velocity is not sufficient for a circular ' ...
    'orbit around the Earth at an altitude of 100km\n']);

% Part 2
h_F9 = 100; % km
v_F9 = 9100; % km/hr

r_F9 = R_earth + h_F9; 
v_F9_kms = v_F9/3600; % km/s

fprintf('Falcon 9 Velocity at 100km Altitude: %.3f km/s\n', v_F9_kms);
fprintf(['The Falcon 9 velocity is not sufficient for a circular ' ...
    'orbit around the Earth at an altitude of 100km\n']);

% Part 3
epsilon_NS = ((v_NS_apogee)^2/2)-(mu_earth/r_NS_apogee); % km^2/s^2
epsilon_F9 = ((v_F9_kms)^2/2)-(mu_earth/r_F9); % km^2/s^2

fprintf('New Shepard specific mechanical energy: %.3e km^2/s^2\n', epsilon_NS);
fprintf('Falcon 9 specific mechanical energy: %.3e km^2/s^2\n', epsilon_F9);
fprintf('The New Shepard specific mechanical energy is greater\n')

% Part 4
r_0_100 = R_earth + 100; % km
r_0_1000 = R_earth + 1000; % km

v_0_100 = sqrt(2*mu_earth/r_0_100); % km/s
v_0_1000 = sqrt(2*mu_earth/r_0_1000); % km/s

fprintf('Escape velocity at an altitude of 100km: %.3f km/s\n', v_0_100);
fprintf('Escape velocity at an altitude of 1000km: %.3f km/s\n', v_0_1000);
fprintf(['As altitude increases, escape velocity decreases because a higher ' ...
    'altitude reduces Earths gravitational influence \n']);

%% Problem 4
disp('')
disp('Problem 4')
r_a_earth = 1.0167*AU; % km
r_p_earth = 0.9833*AU; % km

r_a_67P = 5.0829*AU; % km
r_p_67P = 1.0432*AU; % km

% Part 1
e_earth = (r_a_earth - r_p_earth) / (r_a_earth + r_p_earth);
e_67P = (r_a_67P - r_p_67P) / (r_a_67P + r_p_67P);
fprintf('Eccentricity of Earth orbit: %.3f \n', e_earth);
fprintf('Eccentricity of 67P orbit: %.3f \n', e_67P);

% Part 2
a_earth = r_a_earth / (1+e_earth); % km
a_67P = r_a_67P / (1+e_67P); % km

v_earth_max = sqrt(mu_sun * (2 / r_p_earth - 1 / a_earth)); % km/s
v_earth_min = sqrt(mu_sun * (2 / r_a_earth - 1 / a_earth)); % km/s
v_67P_max = sqrt(mu_sun * (2 / r_p_67P - 1 / a_67P)); % km/s
v_67P_min = sqrt(mu_sun * (2 / r_a_67P - 1 / a_67P)); % km/s
fprintf('Maximum inertial velocity of the Earth: %.3f km/s\n', v_earth_max);
fprintf('Minimum inertial velocity of the Earth: %.3f km/s\n', v_earth_min);
fprintf('Maximum inertial velocity of 67P: %.3f km/s\n', v_67P_max);
fprintf('Minimum inertial velocity of 67P: %.3f km/s\n', v_67P_min);

% Part 3
T_earth = 2*pi*sqrt(a_earth^3/mu_sun)/86400; % sec --> mean solar days
T_67P = 2*pi*sqrt(a_67P^3/mu_sun)/86400; % sec --> mean solar days
fprintf('Orbital period of the Earth: %.3f mean solar days\n', T_earth);
fprintf('Orbital period of 67P: %.3f mean solar days\n', T_67P);