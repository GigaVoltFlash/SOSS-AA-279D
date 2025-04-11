%% Space Mechanics Homework 6
% Tycho Bogdanowitsch
clc; clear;
% Constants
R_earth = 6378; % km
DU_earth = R_earth; % km
AU = 1.496e8; % km
mu_earth = 3.986e5; % km^3/s^2
mu_mars = 4.283e4; % km^3/s^2
mu_sun = 1.327e11; % km^3/s^2
R_sun = 6.963e5; % km

%% Problem 2
disp('Problem 2')
% Part b
e = 0;
RAAN = 180;
d = 86400; % seconds in a day
J_2 = 0.001082; 
num_orbits = 60; % orbits
num_days = 4; % days

T_sidereal_day = 86164.1; % seconds

T_solar_day = 86400; % seconds

a = (mu_earth*(num_days*T_solar_day/(2*pi*num_orbits))^2)^(1/3);

T = 2*pi*sqrt(a^3/mu_earth);

n = 2*pi/T;

year_sec = 365.25*24*3600; % days --> seconds
dRAANdt = deg2rad(360)/year_sec;

i = acosd(dRAANdt*(-2/3)*(1/(n*J_2))*((a*(1-e^2)/R_earth))^2);

e_test = 0.9885;
i_test = 63.4;

dRAANdt_test = (3/2)*n*J_2*(R_earth/(a*(1-e_test^2))^2*cosd(i_test));

fprintf('a: %.3f km\n', a);
fprintf('e: %.3f \n', e);
fprintf('i: %.3f degrees\n', i);
fprintf('RAAN: %.3f degrees\n', RAAN);
fprintf('omega: undefined\n');

% Part c
solar_day = 86400; % seconds
orbits_per_day = (solar_day)/(T);
min_swath = 2*pi*R_earth/orbits_per_day;

fprintf('Min swath: %.3f km\n', min_swath);

%% Problem 3
disp('Problem 3')
J_2_mars = 1.9643e-3; 
R_mars = 3397.2; % km
mu_mars = 4.305e4; % km^3/s^2
omega_mars = 7.088e-5; % rad/s
e = 0;
i = 0;

a_initial = (mu_mars/omega_mars^2)^(1/3);

n_0 = omega_mars;

n_goal = omega_mars;

n_initial = n_0*(1+3*J_2_mars*(R_mars/a_initial)^2);

dn = n_goal-n_initial;

dadn = -(2/3)*(a_initial)/(n_goal);

da = dadn*dn;

a_goal = a_initial + da;

fprintf('Initial altitude: %.3f km\n', a_initial);
fprintf('Delta altitude: %.3f km\n', da);
fprintf('Goal altitude: %.3f km\n', a_goal);