% This function is specifically designed to work with no change in delta da or delta dlambda, and only 
% change delta de.
function [delta_v_vals, delta_v_times] = triplet_burn_tangential(delta_da, delta_de, delta_dlambda, SV1_oe_init, u_SV1_init, t0, tf)
    global J2 R_earth mu_earth
    a = SV1_oe_init(1);
    ecc = SV1_oe_init(2);
    i_SV1_init = SV1_oe_init(3);
    i = deg2rad(i_SV1_init);
    eta = sqrt(1 - ecc^2);
    kappa = 3/4*J2*R_earth^2 * sqrt(mu_earth)/(a^(7/2) * eta^4);
    P = 3*(cos(i))^2 - 1;
    Q = 5*(cos(i))^2 - 1;
    n = sqrt(mu_earth/a^3);
    u0 = deg2rad(u_SV1_init); % Converting initial mean argument of latitude to rad.
    omega_dot = kappa * Q;

    % This result is going to be based on the  
    % Impulsive Maneuvers for Formation Reconfiguration Using Relative Orbital Elements, G.Gaias and S.D'Amico
    % method number 12.
    t1 = -100;
    t2 = -100;
    t3 = -100;
    k = 0;
    while ((t1 < t0) || (t2 < t0)) || (t3 < t0)
        k1 = 1 + k;
        k2 = 2 + k;
        k3 = 3 + k;
        Ubar = atan2(delta_de(2),delta_de(1));
        b = sin(Ubar);
        B = delta_de(2);
        uF = convert_t_to_u(tf);
        c = omega_dot/(n + kappa*(eta*P + Q));
        u1 = 1/(1-c)*(Ubar + k1*pi - c*uF);
        u2 = 1/(1-c)*(Ubar + k2*pi - c*uF);
        u3 = 1/(1-c)*(Ubar + k3*pi - c*uF);
        q = uF - u1 - c/(1 - c)*uF;
        p = uF - u2 - c/(1 - c)*uF;
        l = uF - u3 - c/(1 - c)*uF;
        U1 = (1-c)*u1 + c*uF;
        U2 = (1-c)*u2 + c*uF;
        U3 = (1-c)*u3 + c*uF;
        A = delta_da;
        D = 12*b*(q - l);
        delta_dlambda_dash = (n + kappa*(eta*P + Q))/n*delta_dlambda + 3/2*c/(1 -c)*uF*delta_da;
        if sin(U1) > 0
            g = -1;
        else
            g = 1;
        end
        delta_v_t2 = [0, n*a/4 * (b*A + g*B)/b, 0];
        delta_v_t1 = [0, -n*a*(4*b*delta_dlambda_dash + 3*g*(p-l)*B + 3*b*(p+l)*delta_da)/D, 0];
        delta_v_t3 = [0, n*a*(4*b*delta_dlambda_dash - 3*g*(q-p)*B + 3*b*(q+p)*delta_da)/D, 0];
        
        t1 = convert_u_to_t(u1);
        t2 = convert_u_to_t(u2);
        t3 = convert_u_to_t(u3);
        k = k + 1;
    end
    delta_v_vals = -[delta_v_t1; delta_v_t2; delta_v_t3];
    delta_v_times = [t1, t2, t3];

    if any(delta_v_times > tf)
        fprintf("WARNING! There is a maneuver time that is beyond your control window")
    end

    function t = convert_u_to_t(u)
        t = (u - u0) /(n + kappa*(eta*P + Q)) + t0;
    end
    function u = convert_t_to_u(t)
        u = (n + kappa*(eta*P + Q))*(t - t0) + u0;
    end
end

% function [delta_v_vals, delta_v_times] = triplet_burn_tangential(roe_init, roe_final, SV1_oe_init, u_SV1_init, t0, tf)
%     m = length(t_maneuvers);
%     STM_from_init = calc_STM_for_control(t_0, t_final, SV1_oe_init);
%     no_control_roe_final = STM_from_init * roe_init';
%     total_delta_roe = roe_final' - no_control_roe_final; % Total roe change
%     Matrix_block = zeros(6, 3*m);

%     for i=1:m
%         t_val = t_maneuvers(i);
%         STM = calc_STM_for_control(t_val, t_final, SV1_oe_init);
%         Gamma = calc_Gamma_for_control(t_val, SV1_oe_init, u_SV1_init);
%         Matrix_block(:, 3*(i-1)+1:3*i) = STM * Gamma;
%     end
 
%     delta_v_vals = reshape(pinv(Matrix_block) * total_delta_roe,m,3)';
    
% end