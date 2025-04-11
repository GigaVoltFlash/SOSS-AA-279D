%% Space Mechanics Homework 4
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

%% Problem 2
fprintf(['\nProblem 2\n']);
% Part a
UT1 = [2 1.5 2024];

MJD = Cal2MJD(UT1);

fprintf('The MJD is: %.3f days\n', MJD);

function MJD = Cal2MJD(UT1)
    M = UT1(1);
    D = UT1(2);
    Y = UT1(3);
    if M<=2
        y = Y-1;
        m = M+12;
    else
        y = Y;
        m = M;
    end

    if (y<=1582 && m<=8 && D<= 4)
        B = -2 + (y+4716) / 4 - 1179;
    else
        B = y/400 - y/100 + y/4;
    end
     
    MJD = 365*y - 679004 + floor(B) + floor(30.6001 * (m+1))+D;
end

% Part b
GMST = MJD2GMST(MJD);
GMST_rad_wrapped = wrapTo2Pi(GMST*pi/180); 

fprintf('The GMST is: %.3f rad\n', GMST_rad_wrapped);

function GMST = MJD2GMST(MJD)
    d = MJD - 51544.5;
    GMST = 280.4606 + 360.9856476*d;
end

% Part c
Rot = CRF2TRFRot(GMST_rad_wrapped);
fprintf('The CRF to TRF rotation matrix is: \n');
disp(Rot)

function Rot = CRF2TRFRot(GMST)
    Rot = [cos(GMST) sin(GMST) 0;
        -sin(GMST) cos(GMST) 0;
        0 0 1];
end

%% Problem 3
close all;
fprintf(['\nProblem 3\n']);

% Part b
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

function [x,y] = OE2Perifocal(a,E,e)
    nu = 2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e));
    r = a*(1-e^2)/(1+e*cos(nu));
    x = r*cos(nu);
    y = r*sin(nu);
end

a = 6796.3; % km
RAAN = deg2rad(257.7630); % deg --> rad
e = 0.00023; 
omega = deg2rad(195.9983); % deg --> rad
i = deg2rad(51.6412); % deg --> rad
M = deg2rad(240.3224); % deg --> rad
UT1 = [2 1 2024];

epsilon = 1e-10;

T = 2*pi*sqrt(a^3/mu_earth); % seconds
n = 2*pi/T; % rad/s

delta_t_30_100 = 70*24*3600; % days --> seconds

M_0 = wrapTo2Pi(M + n*delta_t_30_100); % rad

stop_time = 1*24*3600; % days --> seconds

x = out.X;
y = out.Y;

x_start = x(1);
y_start = y(1);

x_end = x(end);
y_end = y(end);

figure(1);
hold on;
plot(x, y,'DisplayName','ISS Orbit');
axis equal;
xlabel('PQW X (km)');
ylabel('PQW Y (km)');
title('ISS Orbit in PQW Frame');
scatter(x_start,y_start,'g','filled','DisplayName', 't = 100 days');
scatter(x_end,y_end,'m','filled','DisplayName', 't = 101 days');
grid on;
legend show;
hold off;

% Part c
function r_ECI = OE2ECI(a,E,e,i,RAAN,omega)
    nu = 2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e));
    r = a*(1-e^2)/(1+e*cos(nu));
    r_PQW = [r * cos(nu); r * sin(nu); 0];

    R_PQR_IJK = [cos(RAAN)*cos(omega)-sin(RAAN)*cos(i)*sin(omega),...
        -cos(RAAN)*sin(omega)-sin(RAAN)*cos(i)*cos(omega),... 
        sin(RAAN)*sin(i);
        sin(RAAN)*cos(omega)+cos(RAAN)*cos(i)*sin(omega),... 
        -sin(RAAN)*sin(omega)+cos(RAAN)*cos(i)*cos(omega),... 
        -cos(RAAN)*sin(i);
        sin(i)*sin(omega), sin(i)*cos(omega), cos(i)];

    r_ECI = R_PQR_IJK*r_PQW;
end

r_ECI = out.R_ECI;

[xE,yE,zE] = ellipsoid(0,0,0,R_earth,R_earth,R_earth,20);

figure(2)
hold on;

surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha', 0.5,'DisplayName', 'Earth');

plot3(r_ECI(1,:), r_ECI(2,:), r_ECI(3,:),'r', 'LineWidth', 1.5,'DisplayName', 'ISS Orbit');

scatter3(r_ECI(1,1), r_ECI(2,1), r_ECI(3,1), 100, 'g', 'filled','DisplayName', 't = 100 days');
scatter3(r_ECI(1,end), r_ECI(2,end), r_ECI(3,end), 100, 'm', 'filled','DisplayName', 't = 101 days');

axis equal;
grid on;
legend show;
xlabel('ECI X (km)');
ylabel('ECI Y (km)');
zlabel('ECI Z (km)');
title('ISS Orbit in ECI Frame');
view(3);
hold off;

% Part d


function r_ECEF = ECEF(a,E,e,i,RAAN,omega,UT1,t)
    M = UT1(1);
    D = UT1(2);
    Y = UT1(3);
    t_days = t/86400;
    D = D+t_days;
    if M<=2
        y = Y-1;
        m = M+12;
    else
        y = Y;
        m = M;
    end

    if (y<=1582 && m<=8 && D<= 4)
        B = -2 + (y+4716) / 4 - 1179;
    else
        B = y/400 - y/100 + y/4;
    end
     
    MJD = 365*y - 679004 + floor(B) + floor(30.6001 * (m+1))+D;
    
    d = MJD - 51544.5;
    GMST = 280.4606 + 360.9856476*d;

    GMST_rad_wrapped = wrapTo2Pi(GMST*pi/180); 

    R_ECI_ECEF = [cos(GMST_rad_wrapped) sin(GMST_rad_wrapped) 0;
        -sin(GMST_rad_wrapped) cos(GMST_rad_wrapped) 0;
        0 0 1];
        
    nu = 2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e));
    r = a*(1-e^2)/(1+e*cos(nu));
    r_PQW = [r * cos(nu); r * sin(nu); 0];

    R_PQR_IJK = [cos(RAAN)*cos(omega)-sin(RAAN)*cos(i)*sin(omega),...
        -cos(RAAN)*sin(omega)-sin(RAAN)*cos(i)*cos(omega),... 
        sin(RAAN)*sin(i);
        sin(RAAN)*cos(omega)+cos(RAAN)*cos(i)*sin(omega),... 
        -sin(RAAN)*sin(omega)+cos(RAAN)*cos(i)*cos(omega),... 
        -cos(RAAN)*sin(i);
        sin(i)*sin(omega), sin(i)*cos(omega), cos(i)];

    r_ECI = R_PQR_IJK*r_PQW;

    r_ECEF = R_ECI_ECEF*r_ECI;
end

r_ECEF = out.R_ECEF;

figure(3)
hold on;

surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha', 0.5,'DisplayName', 'Earth');

plot3(r_ECEF(1,:), r_ECEF(2,:), r_ECEF(3,:),'r', 'LineWidth', 1.5,'DisplayName', 'ISS Orbit');

scatter3(r_ECEF(1,1), r_ECEF(2,1), r_ECEF(3,1), 100, 'g', 'filled','DisplayName', 't = 100 days');
scatter3(r_ECEF(1,end), r_ECEF(2,end), r_ECEF(3,end), 100, 'm', 'filled','DisplayName', 't = 101 days');

axis equal;
grid on;
legend show;
xlabel('ECEF X (km)');
ylabel('ECEF Y (km)');
zlabel('ECEF Z (km)');
title('ISS Orbit in ECEF Frame');
view(3);
hold off;

%% Problem 4
fprintf(['\nProblem 4\n']);
R_earth = 6378.1; % km
e_earth = 0.0818; 

function [lambda,psi] = ECEF2Geocentric(ECEF)
    r_x = ECEF(1);
    r_y = ECEF(2);
    r_z = ECEF(3);

    lambda = atan2d(r_y,r_x);
    psi = asind(r_z/sqrt(r_x^2+r_y^2+r_z^2));
end

r_ECEF = [-2195.7 -4669.6 3761.5];

[lambda,psi] = ECEF2Geocentric(r_ECEF);

fprintf('The Geocentric latitude is: %.3f deg\n', psi);
fprintf('The Geocentric longitude is: %.3f deg\n', lambda);

fprintf('In Geocentric coordinates, the spacecraft is passing over the stadium\n');

function [lambda_prime,psi_prime] = ECEF2Geodetic(ECEF)
    r_x = ECEF(1);
    r_y = ECEF(2);
    r_z = ECEF(3);

    R_earth = 6378.1; % km
    e_earth = 0.0818;     
    epsilon = 0.001; % deg
    
    r_xy = sqrt(r_x^2+r_y^2);
    
    lambda_prime = atan2d(r_y,r_x);

    psi_prime = asind(r_z/sqrt(r_x^2+r_y^2+r_z^2));

    while true
        N = R_earth/(sqrt(1-(e_earth^2)*(sind(psi_prime))^2));
        
        psi_prime_next = atan2d(r_z+N*(e_earth^2)*sind(psi_prime),r_xy);
        delta = psi_prime_next - psi_prime;
        
        if abs(delta)<epsilon
            psi_prime = psi_prime_next;
            break;
        end

        psi_prime = psi_prime_next;
    end
end

[lambda_prime,psi_prime] = ECEF2Geodetic(r_ECEF);

fprintf('The Geodetic latitude is: %.3f deg\n', psi_prime);
fprintf('The Geodetic longitude is: %.3f deg\n', lambda_prime);


arc_dist = 2*pi*R_earth*((psi_prime-36.0909)/360);
fprintf('In Geodetic coordinates, the spacecraft is %.3f km away in the Northward direction\n', arc_dist);

%geoc2geod(psi,R_earth*10e3)