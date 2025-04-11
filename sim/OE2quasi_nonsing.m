function [a, e_x, e_y, i, RAAN, u] = OE2quasi_nonsing(a, e, i, RAAN, w, M)
    % All angles are in degrees
    e_x = e * cosd(w);
    e_y = e * sind(w);
    u = w + M;  % Argument of latitude [deg]
end