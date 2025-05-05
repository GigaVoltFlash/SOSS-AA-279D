function states = evaluate_ya(t, a, e, M_init, K)
    % K is the integration constants that define the YA dynamics
    % t is the total true anomaly values that we want to evaluate HCW over
    % a is the semi-major axis of the chief
    % e is the eccentricity of the chief's orbit
    % M_init is the starting true anomaly of the chief's orbit
    global mu_earth;
    n = sqrt(mu_earth/a^3);

    N = length(t);
    states = zeros(N, 6);
    for i=1:N
        M = deg2rad(M_init) + n*t(i);
        %M = mod(M, 2*pi);
        
        %E = Newton_Raphson(M, e, 1e-5); % in radians
        %f = 2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e)); % in radians
        %f = mod(f,2*pi);
        E = mean2ecc(M,e,1e-10);
        f = ecc2true(E,e);

        eta = sqrt(1-e^2);
        tau = n*t(i)/eta^3; 
        k = 1 + e*cos(f);
        k_prime = -e*sin(f);
        
        A = [a*(eta^2)*eye(3), zeros(3); zeros(3), a*(n/eta)*eye(3)];
        B = compute_B(e, tau, f, k, k_prime);
        total = A * B;
        states(i, :) = (total * K)';
    end
end

function B = compute_B(e, tau, f, k, k_prime)
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
    psi_z_dot_6 = -sin(f);

    B = [psi_x_1 psi_x_2 psi_x_3 0 0 0; ...
        psi_y_1 psi_y_2 psi_y_3 psi_y_4 0 0; ...
        0 0 0 0 psi_z_5 psi_z_6; ...
        psi_x_dot_1 psi_x_dot_2 psi_x_dot_3 0 0 0; ...
        psi_y_dot_1 psi_y_dot_2 psi_y_dot_3 psi_y_dot_4 0 0; ...
        0 0 0 0 psi_z_dot_5 psi_z_dot_6];
end