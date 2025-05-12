function [mean_roe_new, new_time] = ...
    roe_state_space(mean_oe_chief, u_RTN_dep, u_RTN_chief, mean_roe_curr, t, old_time)

J2 = 1.08263e-3;
R_EARTH = 6378136300e-2;       % [m]
MU_EARTH = 0.3986004415e15;    % [m^3/s^2]
new_time = t;
simstep = t - old_time;

% Initialization
mean_roe_new = zeros(6,1);
x_curr = mean_roe_curr;

% Chief orbital elements
a     = mean_oe_chief(1);
e     = mean_oe_chief(2);
inc   = mean_oe_chief(3);
OMEGA = mean_oe_chief(4);
omega = mean_oe_chief(5);
M     = mean_oe_chief(6);

PHI = eye(6);
Q2 = diag(ProcessNoiseCovVec);
R = diag(MeasNoiseCovVec);

if t == 0
    x_apriori = Measurement;
end

% State transition matrix PHI
n = sqrt(MU_EARTH / a^3);
eta = sqrt(1 - e^2);

gamma = 3/4 * J2 * (R_EARTH)^2 * sqrt(MU_EARTH);

E = 1 + eta;
F = 4 + 3*eta;
G = 1 / eta^2;

P2 = 3 * cos(inc)^2 - 1;
Q = 5 * cos(inc)^2 - 1;
S = sin(2*inc);
T = sin(inc)^2;

kappa = gamma / a^(7/2) / eta^4;
dot_omega = kappa * Q;
omega_f = omega + dot_omega * simstep;

% Eccentricity vector rotation
ey0 = e * sin(omega);
ex0 = e * cos(omega);
eyf = e * sin(omega_f);
exf = e * cos(omega_f);

% Populate PHI
PHI(2,1) = -7/2 * kappa * E * P2 * simstep - 3/2 * n * simstep;
PHI(3,1) =  7/2 * kappa * eyf * Q * T;
PHI(4,1) = -7/2 * kappa * exf * Q * T;
PHI(6,1) =  7/2 * kappa * S * simstep;

PHI(2,3) =  kappa * ex0 * F * G * P2 * simstep;
PHI(3,3) =  cos(dot_omega * simstep) - 4 * kappa * ex0 * eyf * G * Q * simstep;
PHI(4,3) =  sin(dot_omega * simstep) + 4 * kappa * ex0 * exf * G * Q * simstep;
PHI(6,3) = -4 * kappa * ex0 * G * S * simstep;

PHI(2,4) =  kappa * ey0 * F * G * P2 * simstep;
PHI(3,4) = -sin(dot_omega * simstep) - 4 * kappa * ey0 * eyf * G * Q * simstep;
PHI(4,4) =  cos(dot_omega * simstep) + 4 * kappa * ey0 * exf * G * Q * simstep;
PHI(6,4) = -4 * kappa * ey0 * G * S * simstep;

PHI(2,5) = -kappa * F * S * simstep;
PHI(3,5) =  5 * kappa * eyf * S * simstep;
PHI(4,5) = -5 * kappa * exf * S * simstep;
PHI(6,5) =  2 * kappa * T * simstep;

% Control input matrix B
E = meantoecc(M, e);
f = 0;
theta = ecctotrue(E, e);
eta = f + omega;

e_x = (ex0 + exf)/2;
e_y = (ey0 + eyf)/2;

B = [
    2 / eta * e * sin(f),      2 / eta * (1 + e * cos(f));
   -2 * eta^2 / (1 + e * cos(f)), 0;
    eta * sin(theta) * ((2 + e * cos(f)) * cos(theta) + e_x) / (1 + e * cos(f)), ...
        eta * e_y * sin(theta) / tan(inc);
   -eta * cos(theta) * ((2 + e * cos(f)) * sin(theta) + e_y) / (1 + e * cos(f)), ...
       -eta * e_x * sin(theta) / cos(inc);
    0, 0;
    0, 0;
];

B = (1 / n) * B; % add /a i f s tate i s given dimensionless ! !

x = PHI * x_curr + B * (u_RTN_dep - u_RTN_chief) * simstep;
mean_roe_new = x(:);

% Measurement matrix
%H = eye(6);

% Kalman Time Update

%P = PHI * P_apriori * PHI' + Q2;

% Measurement Update
%z = H * x;
%if mod(t, 10) < 2
%    K = P * H' / (H * P * H' + R);
%    x = x + K * (Measurement - z);
%    P_updated = (eye(6) - K * H) * P * (eye(6) - K * H)' + K * R * K';
%else
    % P_updated = P; % Optional if no update
%end

%y = H * x;


end
