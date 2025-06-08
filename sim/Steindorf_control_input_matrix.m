function B = Steindorf_control_input_matrix(SV1_OE_state)
% STEINDORF_CONTROL_INPUT_MATRIX Constructs the Steindorf control input matrix (QNS form)
%
% Inputs:
%   SV1_OE_state - 1x6 vector: [a [km], ex, ey, i [deg], RAAN [deg], u [deg]]
%
% Output:
%   B - 6x3 control input matrix

    % Unpack the state vector
    a     = SV1_OE_state(1); % km
    ex    = SV1_OE_state(2);
    ey    = SV1_OE_state(3);
    i_deg = SV1_OE_state(4); % deg
    RAAN_deg = SV1_OE_state(5); % deg (not used here)
    u_deg = SV1_OE_state(6); % deg

    % Convert angles to radians
    i     = deg2rad(i_deg);
    u     = deg2rad(u_deg);

    % Calculate true anomaly from arg of latitude, f = u - omega
    f = u - atan2(ey,ex);
    
    % Compute derived quantities
    e     = sqrt(ex^2 + ey^2);
    eta   = sqrt(1 - e^2);
    one_plus_ecosf = 1 + e * cos(f); 
    
    % Mean motion
    global mu_earth;
    %mu_earth = 398600.4418; % km^3/s^2
    n = sqrt(mu_earth / a^3);   % mean motion [rad/s]

    % Build the matrix
    B = zeros(6,3);

    % Row 1
    B(1,1) = (2 / eta) * e * sin(f);
    B(1,2) = (2 / eta) * one_plus_ecosf;
    B(1,3) = 0;

    % Row 2
    B(2,1) = - (2 * eta^2) / one_plus_ecosf;
    B(2,2) = 0;
    B(2,3) = 0;

    % Row 3
    B(3,1) = eta * sin(u);
    B(3,2) = eta * ( (2 + e * cos(f)) * cos(u) + ex ) / one_plus_ecosf;
    B(3,3) = (eta * ey / tan(i)) * (sin(u) / one_plus_ecosf);

    % Row 4
    B(4,1) = - eta * cos(u);
    B(4,2) = eta * ( (2 + e * cos(f)) * sin(u) + ey ) / one_plus_ecosf;
    B(4,3) = - (eta * ex / tan(i)) * (sin(u) / one_plus_ecosf);

    % Row 5
    B(5,1) = 0;
    B(5,2) = 0;
    B(5,3) = eta * cos(u) / one_plus_ecosf;

    % Row 6
    B(6,1) = 0;
    B(6,2) = 0;
    B(6,3) = eta * sin(u) / one_plus_ecosf;

    % Scale by 1/(a * n)
    B = B / (a * n);

end
