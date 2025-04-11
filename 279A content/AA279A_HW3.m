%% Space Mechanics Homework 3
% Tycho Bogdanowitsch
clc; clear;
% Constants
R_earth = 6371; % km
DU_earth = R_earth; % km
AU = 1.496e8; % km
mu_earth = 3.986e5; % km^3/s^2
mu_mars = 4.283e4; % km^3/s^2
mu_sun = 1.327e11; % km^3/s^2
R_sun = 6.963e5; % km

%% Problem 1
close all;
e = 0.89;
r_a = 3800; % kmco
nu_0 = deg2rad(110); % rad

nu = deg2rad(180); % rad
a = r_a/(1+e); % km
n = sqrt(mu_mars/a^3);the re

E_0 = acos((e+cos(nu_0))/(1+e*cos(nu_0)));
E = acos((e+cos(nu))/(1+e*cos(nu)));

t = (1/n)*((E-e*sin(E))-(E_0-e*sin(E_0))) / 60;

fprintf('The time remaining is: %.3f minutes\n', t);

%% Problem 2
close all;
a = 17.834 * AU; % km
e = 0.967;

% Part i
T = 2*pi*sqrt(a^3/mu_sun); % seconds
T_converted = T / 3.154e7; % seconds --> Earth years
n = 2*pi/T; % rad/s

fprintf('The orbital period of Halleys Comet is: %.3f Earth years\n', T_converted);
fprintf('The mean motion of Halleys Comet is: %.3e rad/s\n', n);

%% Problem 2, Part ii
epsilon = 10e-10; % rad

function E = NewtonRaphson(M,e,epsilon)
    E = M;

    while true
        fE = E - e*sin(E) - M;
        fprimeE = 1 - e*cos(E);
        delta = -fE/fprimeE;
        E_next = E+delta;
        if abs(delta)<epsilon
            break;
        end
    
        E = E_next;
    end
end

%% Problem 2, Part iii
a = 17.834 * AU; % km
e = 0.967;
epsilon = 1e-10;

%% Problem 2, Part iv
close all;
a = 17.834 * AU; % km
e = 0.967;
epsilon = 1e-10;
DU_sun = R_sun;
TU_sun = sqrt(DU_sun^3/mu_sun); % 

T = 2*pi*sqrt(a^3/mu_sun); % seconds
n = 2*pi/T; % rad/s
dt = T/1000;
time = 0:dt:T;

x = out.X;
y = out.Y;

n = length(x);
num_points = 25;

points = round(linspace(1,n,num_points));
x_points = x(points);
y_points = y(points);

%assignin('base', 'dt', dt);
%assignin('base', 'a', a);
%assignin('base', 'e', e);
%assignin('base', 'n', n);
%assignin('base', 'epsilon', epsilon);

%open_system('HalleysCometModel');

%simOut = sim('HalleysCometModel');

%outputData = simOut.get('yout');

%t = outputData.time;
%x_timeseries = outputData.get('x').Values;       % Extract 'x' data
%y_timeseries = outputData.get('y').Values;
%x_data = x_timeseries.Data;
%y_data = y_timeseries.Data;
figure(1);
hold on;
plot(x, y,'DisplayName','Comet Orbit');
scatter(x_points,y_points,'r','filled','DisplayName', '0.04T Orbit Traces');
axis equal;
xlabel('x (meters)');
ylabel('y (meters)');
title('Halley''s Comet Orbit');
grid on;
legend show;
hold off;
% Extract 'y' data
%y = outputData.signals.values;
%plot(t, y);
%xlabel('Time (s)');
%ylabel('Orbital Position');
%title('Halley’s Comet Orbit Simulation');

%% Problem 3
close all;
% Part i
DU_earth = 6378; % km (spherical assumption)
TU_earth = sqrt(DU_earth^3/mu_earth); % 
r = [8;0;0]*DU_earth; % km
v = [-1/sqrt(8); 1/sqrt(8); 0]* DU_earth/TU_earth; % km/s
r_norm = norm(r);
v_norm = norm(v);

epsilon = (v_norm^2/2) - (mu_earth/r_norm);
fprintf('The specific mechanical energy is: %.3f km^2/s^2\n', epsilon);

% Part ii
h = cross(r,v);
h_norm = norm(h);

p_debris = h_norm^2/mu_earth;
r_p_debris = p_debris/2;

e = (cross(v,h)/mu_earth) - r/r_norm;

fprintf('The position of periapsis for the debris is: %.3f km\n', r_p_debris);
fprintf('The threshold position is: %.3f km\n', DU_earth+19134);

% Part iii
nu = 0:1:360; % degrees, true anomaly
e_debris = 1;

r_ECI = p_debris ./(1+e_debris*cosd(nu));
x_parabola = -r_ECI.*sind(nu); % shift to get parabola in right orientation
y_parabola = r_ECI.*cosd(nu);

e_sat = 0;
r_p_sat = 19134+DU_earth; % km

a_sat = r_p_sat/(1-e_sat);
p_sat = a_sat*(1-e_sat^2);
r_ECI = p_sat ./(1+e_sat*cosd(nu));
x_circle = r_ECI.*cosd(nu);
y_circle = r_ECI.*sind(nu);

figure(1)
hold on
scatter(x_circle, y_circle, 10, 'g', 'filled', 'DisplayName', 'Satellite');
scatter(x_parabola, y_parabola, 10, 'm', 'filled', 'DisplayName', 'Debris');
scatter(0, 0, 50, 'k', 'filled', 'DisplayName', 'Central Body');
xlabel('X Position (km)')
ylabel('Y Position (km)')
grid on
legend show;
ylim([-2e5 2e5])
axis equal
hold off

% Part iv
r_p_debris = p_debris/2;
v_debris_p = sqrt(mu_earth * (2 / r_p_debris));

v_sat = sqrt(mu_earth * (2 / r_p_sat - 1 / (2*r_p_sat)));

v_rel = v_debris_p - v_sat;

fprintf('The relative velocity is: %.3f km/s in the -i direction\n', v_rel);

nu = acosd((p_debris/r_norm) -1);
fprintf('The true anomaly of the debris at the initial position is: -%.3f degrees\n', nu);

%% Problem 4
clear;
r_IJK = [3105.4128; -880.8531; 6368.9408]; % km
v_IJK = [-6.7276;-0.5612;3.2023]; % km/s

[a,e,i,RAAN,omega,nu] = ECI2OE(r_IJK,v_IJK);

fprintf('Envisat semi-major axis: %.3f km\n', a);
fprintf('Envisat eccentricity: %.3e. \n', e);
fprintf('Envisat inclination: %.3f degrees\n', i);
fprintf('Envisat RAAN: %.3f degrees\n', RAAN);
fprintf('Envisat argument of perigee: %.3f degrees\n', omega);
fprintf('Envisat true anomaly: %.3f degrees\n', nu);

function [a,e,i,RAAN,omega,nu] = ECI2OE(r_IJK,v_IJK)
    mu_earth = 3.986e5; % km^3/s^2
    r_norm = norm(r_IJK);
    r_i = r_IJK(1);
    r_j = r_IJK(2);
    r_k = r_IJK(3);
    v_norm = norm(v_IJK);
    h = cross(r_IJK,v_IJK);
    h_norm = norm(h);
    W = h/norm(h);
    W_i = W(1);
    W_j = W(2);
    W_k = W(3);
    i = atan2(sqrt(W_i^2 + W_j^2),W_k);
    RAAN = atan2(W_i,-W_j);
    p = h_norm^2/mu_earth;
    a = ((2/r_norm)-(v_norm^2/mu_earth))^(-1);
    n = sqrt(mu_earth/a^3);
    e = sqrt(1-p/a);
    E = atan2(dot(r_IJK,v_IJK)/(a^2*n),(1-r_norm/a));
    M = E - e*sin(E);
    nu = 2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e));
    u = atan2(r_k/sin(i),r_i*cos(RAAN)+r_j*sin(RAAN));
    omega = u - nu;
    i = rad2deg(i);
    RAAN = rad2deg(RAAN);
    omega = rad2deg(omega);
    nu = rad2deg(nu);
    if nu<0
        nu = 360 + nu;
    end
end