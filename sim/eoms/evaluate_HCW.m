function states = evaluate_HCW(t, a, K)
    % K is the integration constants that define the HCW dynamics
    % t is the total time values that we want to evaluate HCW over
    % a is the semi-major axis of the chief
    global mu_earth;
    n = sqrt(mu_earth/a^3);
    A = [a*eye(3), zeros(3); zeros(3), a*n*eye(3)];
    N = length(t);
    states = zeros(N, 6);
    for i=1:N
        B = compute_B(n, t(i));
        total = A * B;
        states(i, :) = (total * K)';
    end
end

function B = compute_B(n, t)
    B = [1, sin(n*t), cos(n*t), 0, 0, 0; ...
        -3/2*n*t, 2*cos(n*t), -2*sin(n*t), 1, 0, 0; ...
        0, 0, 0, 0, sin(n*t), cos(n*t); ...
        0, cos(n*t), -sin(n*t), 0, 0, 0; ...
        -3/2, -2*sin(n*t), -2*cos(n*t), 0, 0, 0; ...
        0, 0, 0, 0, cos(n*t), -sin(n*t)];
end