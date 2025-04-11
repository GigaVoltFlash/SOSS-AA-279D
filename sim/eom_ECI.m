function [statedot] = eom_ECI(t, state, withJ2)
    mu = 3.986e5; % km^3/s^2
    J2 = 1.08263e-3;
    R_earth = 6378.13; % km
    
    % state vector:[rx ry rz vx vy vz]'
    statedot = zeros(6, 1);
    statedot(1:3) = state(4:6); % rdot = velocity
    
    r = state(1:3);
    r_norm = norm(r);
    accel_g = -(mu/r_norm^3)*r;
    a = accel_g;
    
    if withJ2
        r_x = state(1);
        r_y = state(2);
        r_z = state(3);
    
        accel_J2 = (3*J2*mu*R_earth^2/(2*r_norm^5))*...
        [((5*r_z^2/r_norm^2)-1)*r_x;
            ((5*r_z^2/r_norm^2)-1)*r_y;
            r_z*((5*r_z^2/r_norm^2)-3)];

        a = accel_g+accel_J2;
    end

    statedot(4:6) = a;
end 