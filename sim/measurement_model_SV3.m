function y_pred = measurement_model_SV3(x_SV3_pred,SV1_OE)
    % Inputs:
    % x_SV3_pred - predicted scaled ROE of SV3
    % SV1_OE - absolute qns OE of SV1 [a, e_x, e_y, i, RAAN, u] in degrees

    d_a = x_SV3_pred(1);
    d_lambda = x_SV3_pred(2);
    d_e_x = x_SV3_pred(3);
    d_e_y = x_SV3_pred(4);
    d_i_x = x_SV3_pred(5);
    d_i_y = x_SV3_pred(6);

    a = SV1_OE(1);
    ex = SV1_OE(2);
    ey = SV1_OE(3);
    i = SV1_OE(4);
    RAAN = SV1_OE(5);
    u = SV1_OE(6);

    % [a,e,i,RAAN,w,nu, M] = quasi_nonsing2OE(a, ex, ey, i, RAAN, u)
    [a,e,i,RAAN,omega,nu,M] = quasi_nonsing2OE(a, ex, ey, i, RAAN, u);

    % [r_ECI,v_ECI] = ROE2ECI(a_o,e_x_o,e_y_o,i_o,RAAN_o,u_o, ...
    % d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y)
    [r_SV3_ECI,~] = ROE2ECI(a,ex,ey,i,RAAN,u, ...
    d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y);

    chief_oe = [a,e,i,RAAN,omega,M];
    %   chief_oe  = [a, e, i, RAAN, w, M] in degrees
    % [r_ECI,v_ECI] = OE2ECI(a,e,i,RAAN,omega,nu)
    [r_SV1_ECI,v_SV1_ECI] = OE2ECI(a,e,i,RAAN,omega,nu);

    %   roe       = [da; d_lambda; dex; dey; dix; diy] relative orbital elements
    %   chief_oe  = [a, e, i, RAAN, w, M] in degrees
    % Output:
    %   rtn       = 1x3 relative position in RTN frame [km]
    r_SV3_RTN = ROE2RTN(x_SV3_pred, chief_oe, r_SV1_ECI, v_SV1_ECI);

    y_pred = [r_SV3_RTN,r_SV3_ECI']';
end