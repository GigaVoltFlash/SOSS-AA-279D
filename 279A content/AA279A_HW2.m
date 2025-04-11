%% Space Mechanics Homework 2
% Tycho Bogdanowitsch
clc; clear;
% Constants
R_earth = 6371; % km
DU_earth = R_earth; % km
mu_earth = 3.986e5; % km^3/s^2

%% Problem 2
close all;
r_p = 3.0*DU_earth; % km

% Part i (ellipse)
nu = 0:1:360; % degrees, true anomaly
e = 0.7;

a = r_p/(1-e);
p = a*(1-e^2);
r_ECI = p ./(1+e*cosd(nu));
x_ellipse = r_ECI.*cosd(nu);
y_ellipse = r_ECI.*sind(nu);

figure(1)
hold on
scatter(x_ellipse,y_ellipse,10,'blue','filled');

plot([0, r_p], [0, 0], 'r-', 'LineWidth', 2); 
text(r_p/2, 0.1, sprintf('r_p = %.2f km ', r_p), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top','FontSize', 10, 'Color', 'r');

plot([-a*e, -(a+a*e)], [0, 0], 'g--', 'LineWidth', 2); 
text(-(a+a*e)/1.5, -0.1, sprintf('a = %.2f km ', a), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 10, 'Color', 'g');

plot([0, 0], [0, p], 'm-.', 'LineWidth', 2); 
text(0.1, p/2, sprintf('p = %.2f km ', p), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top','FontSize', 10, 'Color', 'm');

plot(0, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'g'); 
text(0, 0, 'Central Body', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 10);
title('Shape of orbit with eccentricity of 0.7')
xlabel('X Position (km)')
ylabel('Y Position (km)')
grid on
axis equal
hold off

% Part ii (circle)
nu = 0:1:360; % degrees, true anomaly
e = 0;

a = r_p/(1-e);
p = a*(1-e^2);
r_ECI = p ./(1+e*cosd(nu));
x_circle = r_ECI.*cosd(nu);
y_circle = r_ECI.*sind(nu);

figure(2)
hold on
scatter(x_circle,y_circle,10,'blue','filled');

plot([0, r_p], [0, 0], 'r-', 'LineWidth', 2); 
text(r_p/2, 0.1, sprintf('r_p = %.2f km ', r_p), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top','FontSize', 10, 'Color', 'r');

plot([-a*e, -(a+a*e)], [0, 0], 'g--', 'LineWidth', 2); 
text(-(a+a*e)/1.5, -0.1, sprintf('a = %.2f km ', a), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 10, 'Color', 'g');

plot([0, 0], [0, p], 'm-.', 'LineWidth', 2); 
text(0.1, p/2, sprintf('p = %.2f km ', p), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','FontSize', 10, 'Color', 'm');

plot(0, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'g'); 
text(0, 0, 'Central Body', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 10);
title('Shape of orbit with eccentricity of 0.0')
xlabel('X Position (km)')
ylabel('Y Position (km)')
grid on
axis equal
hold off

% Part iii (hyperbola)
e = 1.6;

delta = rad2deg(2*asin(1/e));
valid_high = 90 + delta/2;
valid_low = 360 - 90 - delta/2;

nu_1 = 0:1:valid_high-1;
nu_2 = valid_low+1:1:360;
nu = [nu_1 nu_2]; % degrees, true anomaly

a = r_p/(1-e);
p = a*(1-e^2);
r_ECI = p ./(1+e*cosd(nu));
x_hyperbola = r_ECI.*cosd(nu);
y_hyperbola = r_ECI.*sind(nu);

x_asymptote = linspace(-1.5e5, 1.5e5, 1000);
slope_asymptote = sqrt(e^2 - 1); 

y_asymptote_positive = slope_asymptote * x_asymptote - (slope_asymptote*(r_p-a));
y_asymptote_negative = -slope_asymptote * x_asymptote + (slope_asymptote*(r_p-a));

figure(3)
hold on
scatter(x_hyperbola,y_hyperbola,10,'blue','filled');

plot(x_asymptote, y_asymptote_positive, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Asymptote (+)');
plot(x_asymptote, y_asymptote_negative, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Asymptote (-)');

plot([0, r_p], [0, 0], 'r-', 'LineWidth', 2); 
text(r_p/2, 0.1, sprintf('r_p = %.2f km ', r_p), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top','FontSize', 10, 'Color', 'r');

plot([r_p,r_p+(-a)], [0, 0], 'g--', 'LineWidth', 2); 
text(-(a+a*e)/1.5, -0.1, sprintf('a = %.2f km ', -a), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 10, 'Color', 'g');

plot([0, 0], [0, p], 'm-.', 'LineWidth', 2); 
text(0.1, p/2, sprintf('p = %.2f km ', p), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top','FontSize', 10, 'Color', 'm');

plot(0, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'g'); 
text(0, 0, 'Central Body', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 10);
title('Shape of orbit with eccentricity of 1.6')
xlabel('X Position (km)')
ylabel('Y Position (km)')
axis equal
xlim([-2e5 2e5])
grid on

hold off

% Part iv (parabola)
nu_1 = 0:1:179;
nu_2 = 181:1:360;
nu = [nu_1 nu_2]; % degrees, true anomaly
nu = 0:1:360;
e = 1;

a = r_p/(1-e);
p = 2*r_p;
r_ECI = p ./(1+e*cosd(nu));
x_parabola = r_ECI.*cosd(nu);
y_parabola = r_ECI.*sind(nu);

figure(4)
hold on
scatter(x_parabola,y_parabola,10,'blue','filled');

plot([0, r_p], [0, 0], 'r-', 'LineWidth', 2); 
text(r_p/2, 0.1, sprintf('r_p = %.2f km ', r_p), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top','FontSize', 10, 'Color', 'r');

%plot([-a*e, -(a+a*e)], [0, 0], 'g--', 'LineWidth', 2); 
%text(-(a+a*e)/1.5, -0.1, sprintf('a = %.2f km ', a), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 10, 'Color', 'g');

plot([0, 0], [0, p], 'm-.', 'LineWidth', 2); 
text(0.1, p/2, sprintf('p = %.2f km ', p), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','FontSize', 10, 'Color', 'm');

plot(0, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'g'); 
text(0, 0, 'Central Body', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 10);
title('Shape of orbit with eccentricity of 1.0')
xlabel('X Position (km)')
ylabel('Y Position (km)')
xlim([-3e5 3e5])
grid on
axis equal
hold off

figure(5);
hold on;
scatter(x_ellipse, y_ellipse, 10, 'b', 'filled', 'DisplayName', 'Ellipse: e = 0.7');
scatter(x_circle, y_circle, 10, 'g', 'filled', 'DisplayName', 'Circle: e = 0');
scatter(x_hyperbola, y_hyperbola, 10, 'r', 'filled', 'DisplayName', 'Hyperbola: e = 1.6');
scatter(x_parabola, y_parabola, 10, 'm', 'filled', 'DisplayName', 'Parabola: e = 1.0');
scatter(0, 0, 50, 'k', 'filled', 'DisplayName', 'Central Body');
xlabel('X Position (km)');
ylabel('Y Position (km)');
title('Overlay of All Orbit Shapes');
legend show;
axis equal;
xlim([-3e5 3e5])
grid on;
hold off;

%% Problem 3
r_ECI = [-8050 1900 2600]; % km
v_ECi = [-7 -2.8 -4]; % km/s
r_mag = norm(r_ECI); % km
v_mag = norm(v_ECi); % km/s

% Part i
epsilon = (v_mag^2/2) - (mu_earth/r_mag); % km^2/s^2
h = cross(r_ECI,v_ECi); % km^2/s
h_mag = norm(h); % km^2/s
e = sqrt(1 + (2*h_mag^2*epsilon/mu_earth^2));

fprintf('The Eccentricity is: %.3f \n', e);

% Part ii
p = h_mag^2/mu_earth; % km
nu = rad2deg(acos((1/e)*((p/r_mag) - 1))); % deg

fprintf('The True Anomaly is: %.3f deg \n', nu);

% Part iii
fpa = rad2deg(acos(h_mag/(r_mag*v_mag))); % deg

fprintf('The Flight Path Angle is: %.3f deg \n', fpa);

% Part iv
a = -mu_earth/(2*epsilon); % km
E = rad2deg(acos((1/e)*(1-(r_mag/a)))); % deg

fprintf('The Eccentric Anomaly is: %.3f deg \n', E);

% Part v
E_rad = deg2rad(E);
M = rad2deg(E_rad - e*sin(E_rad)); % deg

fprintf('The Mean Anomaly is: %.3f deg \n', M);

%% Problem 4
clear;
a = 6897; % km
e = 0.0002422;
i = 28.4666; % deg
RAAN = 205.7537; % deg
omega = 318.2020; % deg
nu = 41.8386; % deg

[r_result,v_result] = OE2ECI(a,e,i,RAAN,omega,nu);

disp('Hubble Position Vector (km):')
disp(r_result)
disp('Hubble Velocity Vector (km/s):')
disp(v_result)


function [r_ECI,v_ECI] = OE2ECI(a,e,i,RAAN,omega,nu)
    mu_earth = 3.986e5; % km^3/s^2
    i = deg2rad(i);
    RAAN = deg2rad(RAAN);
    omega = deg2rad(omega);
    nu = deg2rad(nu);
    %u = omega + nu;
    r = a*(1-e^2)/(1+e*cos(nu));
    r_PQW = [r * cos(nu); r * sin(nu); 0];
    %r_IJK = r * [cos(u)*cos(RAAN)+sin(u)*cos(i)*sin(RAAN); 
    %    cos(u)*sin(RAAN)+sin(u)*cos(i)*cos(RAAN);
    %    sin(u)*sin(i)];

    n = sqrt(mu_earth/a^3);
    E = 2*atan2(sqrt((1-e)/(1+e))*tan(nu/2),1);
    R_PQR_IJK = [cos(RAAN)*cos(omega)-sin(RAAN)*cos(i)*sin(omega),...
        -cos(RAAN)*sin(omega)-sin(RAAN)*cos(i)*cos(omega),... 
        sin(RAAN)*sin(i);
        sin(RAAN)*cos(omega)+cos(RAAN)*cos(i)*sin(omega),... 
        -sin(RAAN)*sin(omega)+cos(RAAN)*cos(i)*cos(omega),... 
        -cos(RAAN)*sin(i);
        sin(i)*sin(omega), sin(i)*cos(omega), cos(i)];
    v_PQW = (a*n/(1-e*cos(E)))*[-sin(E); sqrt(1-e^2)*cos(E); 0];
    %v_IJK = R_PQR_IJK*v_PQW;

    r_ECI = R_PQR_IJK*r_PQW;
    v_ECI = R_PQR_IJK*v_PQW;

end