function [a,e,i,RAAN,w,nu, M_deg] = quasi_nonsing2OE(a, ex, ey, i, RAAN, u) % all values computed and output in degrees unless specified
    w = rad2deg(atan2(ey, ex));
    e = ex/cosd(w);
    M_deg = u - w;
    nu_rad = mean2true(deg2rad(M_deg),e);
    %E = Newton_Raphson(deg2rad(M), e, 1e-5); % in radians
    nu = rad2deg(nu_rad);
end
