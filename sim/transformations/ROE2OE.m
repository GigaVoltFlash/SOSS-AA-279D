% Takes in  chief (observer) orbit elements and ROE and converts to OE of
% deputy

function [a_t, e_t, i_t, RAAN_t, w_t, nu_t] = ROE2OE(a_o,e_o,i_o,RAAN_o,omega_o,nu_o, ...
    d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y)
    
    M_rad = deg2rad(nu_o);
    M_o = rad2deg(M_rad);

    [a_o,e_x_o,e_y_o,i_o,RAAN_o,u_o] = OE2quasi_nonsing(a_o, e_o, i_o, RAAN_o, omega_o, M_o);

    [a_t,e_x_t,e_y_t,i_t,RAAN_t,u_t] = ...
    ROE2quasi_nonsing(a_o,e_x_o,e_y_o,i_o,RAAN_o,u_o, ...
    d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y);
    
    [a_t, e_t, i_t, RAAN_t, w_t, nu_t, M_t] = ...
    quasi_nonsing2OE(a_t, e_x_t, e_y_t, i_t, RAAN_t, u_t);
end