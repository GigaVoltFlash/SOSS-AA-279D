function position = hcw_equations_with_ic(t, a, state0)
    global mu_earth;
    n = sqrt(mu_earth/a^3);
    N = length(t);
    position = zeros(N, 3);
    x0 = state0(1);
    y0 = state0(2);
    z0 = state0(3);
    xdot0 = state0(4);
    ydot0 = state0(5);
    zdot0 = state0(6);

    position(:, 1) = 4*x0 + 2*ydot0/n + xdot0/n*sin(n*t) - (3*x0 + 2*ydot0/n)*cos(n*t);
    position(:, 2) = -(6*n*x0 + 3*ydot0)*t + (y0 - 2*xdot0/n) + (6*x0 + 4*ydot0/n)*sin(n*t) + 2*xdot0/n * cos(n*t);
    position(:, 3) = zdot0/n*sin(n*t) + z0*cos(n*t);
end