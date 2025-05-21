function [statedot] = secular_J2(t, state)
    global mu_earth;
    global J2;
    global R_earth;
    
    % state vector:[a e_x e_y i RAAN u]'
    statedot = zeros(6, 1);

    a = state(1);
    e_x = state(2);
    e_y = state(3);
    i = deg2rad(state(4));

    n = sqrt(mu_earth/a^3);

    statedot(2) = (-3/4)*(n*J2)*((R_earth/(a*(1-(e_x^2+e_y^2))))^2)*e_y*(5*cos(i)^2-1); % de_x/dt
    statedot(3) = (3/4)*(n*J2)*((R_earth/(a*(1-(e_x^2+e_y^2))))^2)*e_x*(5*cos(i)^2-1); % de_y/dt
    statedot(5) = rad2deg((-3/2)*(n*J2)*((R_earth/(a*(1-(e_x^2+e_y^2))))^2)*cos(i)); % dRAAN/dt
    statedot(6) = rad2deg((3/4)*(n*J2)*((R_earth/(a*(1-(e_x^2+e_y^2))))^2)*(sqrt(1-(e_x^2+e_y^2))* ...
        (3*cos(i)^2-1)+(5*cos(i)^2-1))); % du/dt
end 