function rtn = ROE2RTN(roe, chief_oe)
    % Converts relative orbital elements to RTN position using full ROE-to-ECI-to-RTN pipeline
    % Inputs:
    %   roe       = [da; d_lambda; dex; dey; dix; diy] relative orbital elements
    %   a         = chief semi-major axis [m]
    %   n_mean    = mean motion of chief [rad/s]
    %   chief_oe  = [a, e, i, RAAN, w, M] in degrees
    % Output:
    %   rtn       = 1x3 relative position in RTN frame [m]

    % Parse chief orbital elements
    a_o = chief_oe(1);
    e_o = chief_oe(2);
    i_o = chief_oe(3);
    RAAN_o = chief_oe(4);
    w_o = chief_oe(5);
    M_o = chief_oe(6);

    % Convert to quasi-nonsingular elements
    e_x_o = e_o * cosd(w_o);
    e_y_o = e_o * sind(w_o);
    u_o = M_o + w_o;  % mean argument of latitude (deg)

    % Unpack relative elements
    d_a     = roe(1);
    d_lambda= roe(2);
    d_e_x   = roe(3);
    d_e_y   = roe(4);
    d_i_x   = roe(5);
    d_i_y   = roe(6);

    % Get ECI positions
    [r_deputy, v_deputy] = ROE2ECI(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o, ...
                                   d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y);
    [a_c, e_c, i_c, RAAN_c, w_c, nu_c, ~] = quasi_nonsing2OE(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o);
    [r_chief, v_chief] = OE2ECI(a_c, e_c, i_c, RAAN_c, w_c, nu_c);

    % Compute RTN position
    [rho_RTN, ~] = ECI2RTN_rel(r_chief, v_chief, r_deputy, v_deputy);
    rtn = rho_RTN;
end
