function [a,e,i,RAAN,w,nu, M] = quasi_nonsing2OE(a, ex, ey, i, RAAN, u) % all values computed and output in degrees unless specified
    w = rad2deg(atan2(ey, ex));
    e = ex/cosd(w);
    M = u - w;
    E = Newton_Raphson(deg2rad(M), e, 1e-5); % in radians
    nu = rad2deg(2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e)));
end
