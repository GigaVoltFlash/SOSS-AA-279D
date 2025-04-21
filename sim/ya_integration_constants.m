function [K] = ya_integration_constants(initial_RTN, a, e, f)
    % Initial RTN includes the x, y, z, xdot, ydot, zdot that is also used
    % in relative orbital motion calculation
    % a is the semi-major axis of the chief
    global mu_earth;
    f = deg2rad(f);
    n = sqrt(mu_earth/a^3);
    eta = sqrt(1-e^2);
    tau = 0; % initial condition, t = 0
    k = 1 + e*cos(f);
    k_prime = -e*sin(f);
    psi_x_1 = (1/k) + (3/2)*(k_prime*tau);
    psi_x_2 = sin(f);
    psi_x_3 = cos(f);
    psi_y_1 = (-3/2)*k*tau;
    psi_y_2 = (1 + (1/k))*cos(f);
    psi_y_3 = -(1 + (1/k))*sin(f);
    psi_y_4 = (1/k);
    psi_z_5 = (1/k)*sin(f);
    psi_z_6 = (1/k)*cos(f);
    psi_x_dot_1 = (k_prime/2) - (3/2)*(k^2)*(k-1)*tau;
    psi_x_dot_2 = (k^2)*cos(f);
    psi_x_dot_3 = -(k^2)*sin(f);
    psi_y_dot_1 = -(3/2)*(k+(k^2)*k_prime*tau);
    psi_y_dot_2 = -(k^2 + 1)*sin(f);
    psi_y_dot_3 = -e - (k^2 + 1)*cos(f);
    psi_y_dot_4 = -k_prime;
    psi_z_dot_5 = e + cos(f);
    psi_z_dot_6 = sin(f);

    A = [a*(eta^2)*eye(3), zeros(3); zeros(3), (a*n/eta)*eye(3)];
    B = [psi_x_1 psi_x_2 psi_x_3 0 0 0; ...
        psi_y_1 psi_y_2 psi_y_3 psi_y_4 0 0; ...
        0 0 0 0 psi_z_5 psi_z_6; ...
        psi_x_dot_1 psi_x_dot_2 psi_x_dot_3 0 0 0; ...
        psi_y_dot_1 psi_y_dot_2 psi_y_dot_3 psi_y_dot_4 0 0; ...
        0 0 0 0 psi_z_dot_5 psi_z_dot_6];
    total = A * B;
    K = total\initial_RTN;
end