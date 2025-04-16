function [tspan, x_out] = rk4_eom_ECI(tspan, x0, withJ2)
    % RK4 Fixed-step 4th-order Runge-Kutta integrator
    %
    % Inputs:
    %   tspan  - [t0 tf] integration time range with fixed time-step
    %   x0     - initial condition (column or row vector)
    %   withJ2 - Flag when true turns on J2 perturbations
    %
    % Outputs:
    %   tspan  - vector of time values
    %   x_out  - solution matrix (rows correspond to time steps)
    dt = tspan(2) - tspan(1);
    n_steps = length(tspan);
    n_vars = length(x0);
    x_out = zeros(n_steps, n_vars);
    
    x = x0(:);  % ensure column vector
    x_out(1,:) = x.';
    
    for i = 2:n_steps
        t = tspan(i-1);
        
        k1 = eom_ECI(t, x, withJ2);
        k2 = eom_ECI(t + dt/2, x + dt/2 * k1, withJ2);
        k3 = eom_ECI(t + dt/2, x + dt/2 * k2, withJ2);
        k4 = eom_ECI(t + dt,   x + dt   * k3, withJ2);
        
        x = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        x_out(i,:) = x.';
    end
end
