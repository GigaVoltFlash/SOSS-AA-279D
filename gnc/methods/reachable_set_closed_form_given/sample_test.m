% References: Chernick PhD
% Starter code is NOT self-encompassing. You are expected to implement some
% additional functions. Everything is `impulsive_control.m` should have
% function headers at a minimum though. 

close all; clear; clc;

addpath('utils/');
addpath('mean_osc/');

% Constants
J2 = 0.0010826358191967; % Earth's J2 coefficient
mu = 3.986004415e14; % Earth gravitational parameter [m^3/s^2]
rE = 6.378136300e6; % Earth radius [m]

% Plotting
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesFontSize',18);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Integration
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial chief absolute singular orbital parameters
sma_c  = 6892.927e3;      % semi-major axis [m] %sma_c = 10000e3;
ecc_c  = 1e-4;            % eccentricity component in I
inc_c  = deg2rad(97.44);  % inclination [rad]
raan_c = deg2rad(270);    % RAAN [rad]
aop_c  = deg2rad(0);      % aop [rad]
M_c    = deg2rad(0);      % mean anomaly [rad]

oe_init_c = [sma_c, ecc_c, inc_c, raan_c, aop_c, M_c]; % combine the above into a single vector
oe_init_c_j2 = mean2osc(oe_init_c, 1);

[r_init_c,v_init_c] = singular_oe2rv(mu,oe_init_c_j2); % chief position and velocity in inertial frame, J2-perturbed
rv_init_c = [r_init_c;v_init_c];

% Helpful parameters
T = 2*pi*sqrt(sma_c^3/mu); % [sec]
n = sqrt(mu/sma_c^3); % mean motion

% Initial conditions of the Deputy as QNS ROEs [dsma;dlambda;dex;dey;dix;diy]
aroe_init1 = [0, 100, 50,100,30,200]; % m
roe_init1 = aroe_init1 / sma_c;

% Initial absolute conditions of the deputy
oe_init_d1 = roe2singular_oe(oe_init_c,roe_init1);
oe_init_d1_j2 = mean2osc(oe_init_d1, 1);

[r_init_d,v_init_d] = singular_oe2rv(mu,oe_init_d1_j2); % deputy position and velocity in inertial frame, J2 perturbed
rv_init_d = [r_init_d;v_init_d];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconfiguration Conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desired ROE state 
aroe_des = [0, 100, -100, 100, 50, 100]; % m

% Desired reconfiguration time
dt = 10*T; % 10 orbit periods

% Simulation step size
sim_step = 10; % s

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine optimal maneuver times and magnitudes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the functions you'd like to use for control input matrix and STM
stm = @(chief_oe, t) chernick_J2_stm(chief_oe, t, rE, mu, J2); 
control_input_matrix = @(chief_oe) chernick_control_matrix(chief_oe, mu);

% Get the optimal maneuver plan for the reconfiguration window
[t_maneuvers, manuevers, total_cost] = impulsive_control(oe_init_c', aroe_init1', aroe_des', dt, stm, control_input_matrix, rE, mu, J2);

% Round the time
t_maneuvers = round(t_maneuvers./sim_step).*sim_step;