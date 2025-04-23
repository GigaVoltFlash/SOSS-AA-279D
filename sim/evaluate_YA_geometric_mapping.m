function states = evaluate_YA_geometric_mapping(t, a, ex, ey, i, w, M_init, d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y)
    global mu_earth;
    n = sqrt(mu_earth / a^3);
    e = sqrt(ex^2 + ey^2);
    eta = sqrt(1 - e^2);
    

    ROE_init = [d_a;d_lambda;d_e_x;d_e_y;d_i_x;d_i_y];
    N = length(t);
    states = zeros(N, 6);
    for idx = 1:N
        M = deg2rad(M_init) + n * t(idx);
        %M = mod(M, 2*pi);
        E = Newton_Raphson(M, e, 1e-5);
        f = 2 * atan2(sqrt(1 + e) * tan(E / 2), sqrt(1 - e));
        f = mod(f,2*pi);
        u = w + f;  % True longitude
        %u = f;

        %k = 1 + e*cos(f);
        %k_prime = -e*sin(f);

        k = 1 + ex*cos(u)+ ey*sin(u);
        k_prime = -ex*sin(u) + ey*cos(u);

        %tau = n * t(idx) / eta^3;

        A = [a * (eta^2)*eye(3), zeros(3); zeros(3), a * (n / eta) * eye(3)];
        B = compute_B_full(t(idx), k, k_prime, n, eta, u, i, ex, ey);
        total = A * B;
        states(idx, :) = (total * ROE_init)';
    end
end

function B = compute_B_full(t, k, kp, n, eta, u, i, ex, ey)
    
    % Bx block
    bx1 = 1/k + (3/2)*kp*(n/eta^3)*t;
    bx2 = -kp / eta^3;
    bx3 = (1/eta^3)*(ex*(k - 1)/(1 + eta) - cos(u));
    bx4 = (1/eta^3)*(ey*(k - 1)/(1 + eta) - sin(u));
    bx6 = (kp / eta^3) * cot(i);
    
    % By block
    by1 = -(3/2)*k*(n/eta^3) * t;
    by2 = k / eta^3;
    by3 = (1/eta^2)*((1 + 1/k)*sin(u) + (ey/k) + (k/eta)*(ey/(1 + eta)));
    by4 = -(1/eta^2)*((1 + 1/k)*cos(u) + (ex/k) + (k/eta)*(ex/(1 + eta)));
    by6 = (1/k - (k/eta^3))*cot(i);
    
    % Bz block
    bz5 = (1/k) * sin(u);
    bz6 = -(1/k) * cos(u);

    % Bxdot block
    bx_dot1 = (kp/2) + (3/2)*(k^2)*(1 - k)*(n/eta^3) * t;
    bx_dot2 = (k^2 / eta^3) * (k - 1);
    bx_dot3 = (k^2 / eta^3) * (eta * sin(u) + ey * (k - 1)/(1 + eta));
    bx_dot4 = -(k^2 / eta^3) * (eta * cos(u) + ex * (k - 1)/(1 + eta));
    bx_dot6 = -(k^2 / eta^3) * (k - 1) * cot(i);

    % Bydot block
    by_dot1 = -(3/2)*k*(1 + k * kp * (n / eta^3) * t);
    by_dot2 = (k^2 / eta^3) * kp;
    by_dot3 = (1 + (k^2 / eta^3))*cos(u) + ex*(k/eta^2)*(1 + (k/eta)*(1 - k)/(1 + eta));
    by_dot4 = (1 + (k^2 / eta^3))*sin(u) + ey*(k/eta^2)*(1 + (k/eta)*(1 - k)/(1 + eta));
    by_dot6 = -(1 + (k^2 / eta^3))*kp*cot(i);

    % Bz dot block
    bz_dot5 = cos(u) + ex;
    bz_dot6 = sin(u) + ey;

    B = zeros(6, 6);
    B(1,1) = bx1;    B(1,2) = bx2;    B(1,3) = bx3;    B(1,4) = bx4;    B(1,6) = bx6;
    B(2,1) = by1;    B(2,2) = by2;    B(2,3) = by3;    B(2,4) = by4;    B(2,6) = by6;
    B(3,5) = bz5;    B(3,6) = bz6;
    B(4,1) = bx_dot1;B(4,2) = bx_dot2;B(4,3) = bx_dot3;B(4,4) = bx_dot4;B(4,6) = bx_dot6;
    B(5,1) = by_dot1;B(5,2) = by_dot2;B(5,3) = by_dot3;B(5,4) = by_dot4;B(5,6) = by_dot6;
    B(6,5) = bz_dot5;B(6,6) = bz_dot6;
end
