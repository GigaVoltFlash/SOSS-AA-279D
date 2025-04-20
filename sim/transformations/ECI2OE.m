% Converting ECI to Orbital Elements
function [a,e,i,RAAN,omega,nu] = ECI2OE(r_IJK,v_IJK)
    global mu_earth
    r_norm = norm(r_IJK);
    r_i = r_IJK(1);
    r_j = r_IJK(2);
    r_k = r_IJK(3);
    v_norm = norm(v_IJK);
    h = cross(r_IJK,v_IJK);
    h_norm = norm(h);
    W = h/norm(h);
    W_i = W(1);
    W_j = W(2);
    W_k = W(3);
    i = atan2(sqrt(W_i^2 + W_j^2),W_k);
    RAAN = atan2(W_i,-W_j);
    p = h_norm^2/mu_earth;
    a = ((2/r_norm)-(v_norm^2/mu_earth))^(-1);
    n = sqrt(mu_earth/a^3);
    e = sqrt(1-p/a);
    E = atan2(dot(r_IJK,v_IJK)/(a^2*n),(1-r_norm/a));
    M = E - e*sin(E);
    nu = 2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e));
    u = atan2(r_k/sin(i),r_i*cos(RAAN)+r_j*sin(RAAN));
    omega = u - nu;
    i = rad2deg(i);
    RAAN = rad2deg(RAAN);
    omega = rad2deg(omega);
    nu = rad2deg(nu);
    if nu<0
        nu = 360 + nu;
    end
    if omega<-180
        omega = 360 + omega;
    end
end