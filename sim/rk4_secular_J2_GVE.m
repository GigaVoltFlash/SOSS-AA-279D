function [tspan, x_out] = rk4_secular_J2_GVE(tspan, x0)
    % RK4 Fixed-step 4th-order Runge-Kutta integrator
    %
    % Inputs:
    %   tspan  - [t0 tf] integration time range with fixed time-step
    %   x0     - initial condition (column or row vector)
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
        
        k1 = secular_J2(t, x);
        k2 = secular_J2(t + dt/2, x + dt/2 * k1);
        k3 = secular_J2(t + dt/2, x + dt/2 * k2);
        k4 = secular_J2(t + dt,   x + dt   * k3);
        
        x = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        x_out(i,:) = x.';
    end
end
