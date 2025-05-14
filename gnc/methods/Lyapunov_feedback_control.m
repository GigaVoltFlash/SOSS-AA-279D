function control_vec = Lyapunov_feedback_control(roe_curr, roe_desired, oe_chief, N, k, n_orbits)
    % Inputs:
    % roe_curr: current ROE unscaled by a_c
    % roe_desired: desired ROE unscaled by a_c
    % oe_chief: current chief absolute OE
    % N:  4?                   % Must be even and > 2
    % k:   1000?              % Large scalar for scaling
    
    % Output:
    % control_vec: 3x1 control output in RTN m/s^2

    % Useful global constants (in km)
    global J2 R_earth mu_earth

    a_c = oe_chief(1); e_c = oe_chief(2); i_c = deg2rad(oe_chief(3));
    omega_c = deg2rad(oe_chief(5)); nu_c = deg2rad(oe_chief(6));

    a_c_meters = a_c*1e3;

    % using reduced model (remove d_lambda)
    roe_curr_reduced = [roe_curr(1),roe_curr(3),roe_curr(4),roe_curr(5),roe_curr(6)];
    roe_desired_reduced = [roe_desired(1),roe_desired(3),roe_desired(4),roe_desired(5),roe_desired(6)];

    delta_roe = (roe_curr_reduced-roe_desired_reduced); %./a_c; % Unscale ROE if needed
    delta_dlambda = roe_curr(2) - roe_desired(2);
    % Find the A Matrix (Steindorf Eq A.2)
    eta = sqrt(1 - e_c^2);
    gamma = (3/4) * J2 * R_earth^2 * sqrt(mu_earth); % 
    kappa = gamma / (a_c^(7/2) * eta^4);
    
    ex = e_c * cos(omega_c);
    ey = e_c * sin(omega_c);
    
    Q = 5 * cos(i_c)^2 - 1;
    S = sin(2 * i_c);
    T = sin(i_c)^2;
    G = 1 / eta^2;    

    A = zeros(5, 5);
    
    A(2,1) = (7/2) * ey * Q;
    A(2,2) = -4 * ex * ey * G * Q;
    A(2,3) = -(1 + 4 * ey^2 * G) * Q;
    A(2,4) = 5 * ey * S;
    
    A(3,1) = -(7/2) * ex * Q;
    A(3,3) = (1 + 4 * ex^2 * G) * Q;
    A(3,4) =  4 * ex * ey * G * Q;
    A(3,5) = -5 * ex * S;
    
    A(5,1) = (7/2) * S;
    A(5,2) = -4 * ex * G * S;
    A(5,3) = -4 * ey * G * S;
    A(5,4) = 2 * T;

    A = kappa * A;
    
    

    % Compote B Matrix
    f = nu_c;
    theta = f + omega_c;
    n = sqrt(mu_earth/ a_c^3);
    
    cosf = cos(f);
    costh = cos(theta);
    sinth = sin(theta);
    tan_i = tan(i_c);
    
    denom = 1 + e_c * cosf;
    
    B = zeros(5, 2);
    
    B(1,1) = (2 / eta) * denom;
    B(2,1) = eta * ((2 + e_c * cosf) * costh + ex) / denom;
    B(2,2) = eta * ((ey / tan_i) + (sinth / denom));
    B(3,1) = eta * ((2 + e_c * cosf) * sinth + ey) / denom;
    B(3,2) = -eta * ((ex / tan_i) - (sinth / denom));
    B(4,1) = 0;
    B(4,2) = eta * costh / denom;
    B(5,1) = 0;
    B(5,2) = eta * sinth / denom;
    
    B = B / (a_c * n); % or in meters?

    % Compute P Matrix
    M = true2mean(nu_c, e_c); % rad

    u = M + omega_c;

    delta_d_e_y = delta_roe(3);
    delta_d_e_x = delta_roe(2);
    delta_d_i_y = delta_roe(5);
    delta_d_i_x = delta_roe(4);

    u_ip = atan2(delta_d_e_y, delta_d_e_x); % Optimal in-plane argument
    u_oop = atan2(delta_d_i_y, delta_d_i_x); % Optimal out-of-plane argument
    
    J = u - u_ip;
    H = u - u_oop;
    
    % Calculations required for the delta lambda manipulation
    v_opt = a_c*n*norm([delta_d_e_x, delta_d_e_y])/(2*eta);
    delta_v_2pi = v_opt/n_orbits;
    delta_d_a_tan = 2/(a_c*n*eta)*(1 + e_c*cos(f))*delta_v_2pi/2;
    delta_d_a_des = 1e3*abs(delta_d_a_tan)/2;
    delta_d_lambda_rate_des = 3/2*n*abs(delta_d_a_des);
    tau = 1e7;
    if delta_dlambda >= 0
        delta_lambda_rate = -min([abs(delta_dlambda)/tau, delta_d_lambda_rate_des]);
    else
        delta_lambda_rate = min([abs(delta_dlambda)/tau, delta_d_lambda_rate_des]);
    end

    
    delta_da_command =  -2/3*delta_lambda_rate/n*1e3;

    if abs(delta_da_command) >= abs(delta_d_a_des)
        P_diag = [0, ...
              0, ...
              0, ...
              cos(H)^N, ...
              cos(H)^N];
    else
        P_diag = [cos(J)^N, ...
              cos(J)^N, ...
              cos(J)^N, ...
              cos(H)^N, ...
              cos(H)^N];
    end

    P = diag(P_diag) / k;
    
    % Changing delta_da for the manipulation
    delta_roe(1) = roe_curr(1) - delta_da_command;

    
    % Compute control input
    u_2 = -pinv(B) * (A * roe_curr_reduced' + P * delta_roe');

    control_vec = [0; u_2(:)];  % 3x1 vector: [0; u_2(1); u_2(2)]
end
