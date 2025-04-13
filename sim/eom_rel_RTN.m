 function [statedot] = eom_rel_RTN(t, state)
   global mu_earth
   % Extract variables from the main full state
   [state0, x, y, z, xdot, ydot, zdot] = state(1:6), state(7), state(8), state(9), state(10), state(11), state(12);
   % Derivative of the ECI co-ordinates of the chief
   statedot(1:6) = eom_ECI(t, state0, False); % No perturbations in the chief's state
   
   % For the chief's current position, get the angular velocity
   theta_dot = norm(cross(state(1:3), state(4:6))) ./ norm(state(1:3)).^2;  % each row is ω_eci (in ECI)

   % Get the angular acceleration using r0_dot from the RTN conversion
   % (maybe there is a better way for this?)
   r0 = norm(state(1:3));
   [r_RTN, v_RTN] = ECI2RTN(state(1:3), state(4:6));
   r0_dot = v_RTN(1);
   theta_ddot = -2*r0_dot*theta_dot/r0;
   
   statedot(10) = 2*theta_dot*ydot + theta_ddot*y * theta_dot^2*x -mu_earth*(r0 + x)/((r0 + x)^2 + y^2 + z^2)^(3/2) + mu_earth/r0^2;
   statedot(11) = -2*theta_dot*xdot - theta_ddot*x + theta_dot^2*y -mu_earth*y/((r0 + x)^2 + y^2 + z^2)^(3/2);
   statedot(12) = -mu_earth*z/((r0 + x)^2 + y^2 + z^2)^(3/2);

   statedot(7) = xdot;
   statedot(8) = ydot;
   statedot(9) = zdot;

 end

