function [r_ECI,v_ECI] = OE2ECI(a,e,i,RAAN,omega,nu)
    global mu_earth;
    i = deg2rad(i);
    RAAN = deg2rad(RAAN);
    omega = deg2rad(omega);
    nu = deg2rad(nu);
    r = a*(1-e^2)/(1+e*cos(nu));
    r_PQW = [r * cos(nu); r * sin(nu); 0];

    n = sqrt(mu_earth/a^3);
    E = 2*atan2(sqrt((1-e)/(1+e))*tan(nu/2),1);
    R_PQW_IJK = [cos(RAAN)*cos(omega)-sin(RAAN)*cos(i)*sin(omega),...
        -cos(RAAN)*sin(omega)-sin(RAAN)*cos(i)*cos(omega),... 
        sin(RAAN)*sin(i);
        sin(RAAN)*cos(omega)+cos(RAAN)*cos(i)*sin(omega),... 
        -sin(RAAN)*sin(omega)+cos(RAAN)*cos(i)*cos(omega),... 
        -cos(RAAN)*sin(i);
        sin(i)*sin(omega), sin(i)*cos(omega), cos(i)];
    v_PQW = (a*n/(1-e*cos(E)))*[-sin(E); sqrt(1-e^2)*cos(E); 0];

    r_ECI = R_PQW_IJK*r_PQW;
    v_ECI = R_PQW_IJK*v_PQW;
end
