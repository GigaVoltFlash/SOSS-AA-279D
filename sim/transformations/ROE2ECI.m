% Takes in quasi-nonsing chief (observer) elements and ROE and converts to position
% and velocity of the deputy (target)

function[r_ECI,v_ECI] = ROE2ECI(a_o,e_x_o,e_y_o,i_o,RAAN_o,u_o, ...
    d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y)

    J2_flag = 1;

    

    [a_o, e_o, i_o, RAAN_o, w_o, nu_o, M_o] = ...
    quasi_nonsing2OE(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o);

    mean_elem = [a_o*1e3,e_o,deg2rad(i_o),deg2rad(RAAN_o),deg2rad(w_o),deg2rad(M_o)];

    osc_elem = mean2osc(mean_elem, J2_flag);

    a_o_osc = osc_elem(1)/1e3;
    e_o_osc = osc_elem(2);
    i_o_osc = rad2deg(osc_elem(3));
    RAAN_o_osc = rad2deg(osc_elem(4));
    w_o_osc = rad2deg(osc_elem(5));
    nu_o_osc = rad2deg(mean2true(osc_elem(6),e_o_osc));
    M_o_osc = rad2deg(osc_elem(6));

    [a_o_osc, e_x_o_osc, e_y_o_osc, i_o_osc, RAAN_o_osc, u_o_osc] = ...
        OE2quasi_nonsing(a_o_osc, e_o_osc, i_o_osc, RAAN_o_osc, w_o_osc, M_o_osc);
    
    [a_t,e_x_t,e_y_t,i_t,RAAN_t,u_t] = ...
    ROE2quasi_nonsing(a_o_osc,e_x_o_osc,e_y_o_osc,i_o_osc,RAAN_o_osc,u_o_osc, ...
    d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y);

    % [a_t,e_x_t,e_y_t,i_t,RAAN_t,u_t] = ...
    % ROE2quasi_nonsing(a_o,e_x_o,e_y_o,i_o,RAAN_o,u_o, ...
    % d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y);
    
    [a_t, e_t, i_t, RAAN_t, w_t, nu_t, M_t] = ...
    quasi_nonsing2OE(a_t, e_x_t, e_y_t, i_t, RAAN_t, u_t);

%    a - semi-major axis [m]
%                 e - eccentricity [-]
%                 i - inclination [rad]
%                 O - right ascension of the ascending node [rad]
%                 w - argument of periapsis [rad]
%                 M - mean anomaly [rad]

    % mean_elem = [a_t*1e3,e_t,deg2rad(i_t),deg2rad(RAAN_t),deg2rad(w_t),deg2rad(M_t)];
    % 
    % osc_elem = mean2osc(mean_elem, J2_flag);
    % 
    % a_t_osc = osc_elem(1)/1e3;
    % e_t_osc = osc_elem(2);
    % i_t_osc = rad2deg(osc_elem(3));
    % RAAN_t_osc = rad2deg(osc_elem(4));
    % w_t_osc = rad2deg(osc_elem(5));
    % nu_t_osc = rad2deg(mean2true(osc_elem(6),e_t_osc));
    % 
    % [r_ECI, v_ECI] = OE2ECI(a_t_osc, e_t_osc, i_t_osc, RAAN_t_osc, w_t_osc, nu_t_osc);


    osc_elem = [a_t*1e3,e_t,deg2rad(i_t),deg2rad(RAAN_t),deg2rad(w_t),deg2rad(M_t)];
    
    mean_elem = osc2mean(osc_elem, J2_flag);

    a_t_mean = mean_elem(1)/1e3;
    e_t_mean = mean_elem(2);
    i_t_mean = rad2deg(mean_elem(3));
    RAAN_t_mean = rad2deg(mean_elem(4));
    w_t_mean = rad2deg(mean_elem(5));
    nu_t_mean = rad2deg(mean2true(mean_elem(6),e_t_mean));


    [r_ECI,v_ECI] = OE2ECI(a_t_mean, e_t_mean, i_t_mean, RAAN_t_mean, w_t_mean, nu_t_mean);
end