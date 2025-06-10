function [statedot] = secular_J2_and_drag(t, state)
% SECULAR_J2_DRAG Computes secular rates of QNS elements under J2 + drag
% using provided GVE formulas, labeled separately
%
% Inputs:
%   t     - time [s] (unused here)
%   state - 6x1 vector: [a [km], e_x, e_y, i [deg], RAAN [deg], u [deg]]
%
% Outputs:
%   statedot - 6x1 vector of time derivatives of state

    % --- Globals ---
    global mu_earth;
    global J2;
    global R_earth;
    global mass_sc;
    global area_sc;
    global CD_sc;

    % --- Unpack state ---
    a     = state(1);           % km
    e_x   = state(2);
    e_y   = state(3);
    i_deg = state(4);
    RAAN_deg = state(5);        % not used here
    u_deg = state(6);

    % --- Derived quantities ---
    i = deg2rad(i_deg);         % rad
    u_rad = deg2rad(u_deg);     % rad

    e = sqrt(e_x^2 + e_y^2);
    n = sqrt(mu_earth / a^3);   % mean motion [rad/s]
    p = a * (1 - e^2);

    % --- Recover true anomaly ---
    omega = atan2(e_y, e_x);
    f = u_rad - omega;

    % --- Orbital radius ---
    r = a * (1 - e^2) / (1 + e * cos(f)); % km

    % --- Orbital velocity (non-circular) ---
    v = sqrt( mu_earth * (2/r - 1/a) ); % km/s
    v_mps = v * 1e3; % m/s

    % --- Atmospheric density ---
    alt_km = r - R_earth; % current altitude in km
    rho = atmospheric_density_Vallado(alt_km); % kg/m^3

    % --- Drag force ---
    D = 0.5 * rho * v_mps^2 * CD_sc * area_sc; % N

    % --- Acceleration due to drag ---
    accel_drag = D / mass_sc; % m/s^2
    accel_drag_kmps2 = accel_drag / 1e3; % km/s^2

    % --- Perturbing acceleration components ---
    f_R = 0;
    f_T = -accel_drag_kmps2;
    f_N = 0;

    % ----------------------------
    % --- J2 effects ------------
    % ----------------------------
    statedot_J2 = zeros(6,1);

    % de_x/dt (J2)
    statedot_J2(2) = (-3/4) * n * J2 * (R_earth / p)^2 * e_y * (5*cos(i)^2 - 1);

    % de_y/dt (J2)
    statedot_J2(3) = (3/4) * n * J2 * (R_earth / p)^2 * e_x * (5*cos(i)^2 - 1);

    % dRAAN/dt (J2)
    statedot_J2(5) = rad2deg((-3/2) * n * J2 * (R_earth / p)^2 * cos(i));

    % du/dt (J2 + Keplerian)
    du_J2 = (3/4) * n * J2 * (R_earth / p)^2 * (sqrt(1 - e^2)*(3*cos(i)^2 - 1) + (5*cos(i)^2 - 1));
    du_total_J2 = n + du_J2;
    statedot_J2(6) = rad2deg(du_total_J2); % deg/s

    % Note: J2 does not affect a or i:
    statedot_J2(1) = 0;
    statedot_J2(4) = 0;

    % ----------------------------
    % --- Drag effects -----------
    % ----------------------------
    statedot_drag = zeros(6,1);

    % da/dt (drag)
    statedot_drag(1) = (2 / n) * f_T;

    % de_x/dt (drag)
    statedot_drag(2) = (1 / (n * a)) * (sin(u_rad) * f_R + 2 * cos(u_rad) * f_T);

    % de_y/dt (drag)
    statedot_drag(3) = (1 / (n * a)) * (-cos(u_rad) * f_R + 2 * sin(u_rad) * f_T);

    % di/dt (drag)
    statedot_drag(4) = (1 / (n * a)) * cos(u_rad) * f_N;

    % dRAAN/dt (drag)
    statedot_drag(5) = (1 / (n * a * sin(i))) * sin(u_rad) * f_N;

    % du/dt (drag) → no f_T term, only f_R and f_N → both zero here!
    statedot_drag(6) = rad2deg( (1 / (n * a)) * (-2 * f_R - sin(u_rad) * cot(i) * f_N) );

    % ----------------------------
    % --- Combine total ---------
    % ----------------------------

    statedot = statedot_J2 + statedot_drag;

end