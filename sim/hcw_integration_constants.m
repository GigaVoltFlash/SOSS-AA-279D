function [K] = hcw_integration_constants(initial_RTN, a)
    % Initial RTN includes the x, y, z, xdot, ydot, zdot that is also used
    % in relative orbital motion calculation
    % a is the semi-major axis of the chief
    x     = initial_RTN(1);
    y     = initial_RTN(2);
    z     = initial_RTN(3);
    xdot  = initial_RTN(4);
    ydot  = initial_RTN(5);
    zdot  = initial_RTN(6); 

    % a = block()
end