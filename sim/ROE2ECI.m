% Takes in quasi-nonsing chief (observer) elements and ROE and converts to position
% and velocity of the deputy (target)

function[r_ECI,v_ECI] = ROE2ECI(a_o,e_x_o,e_y_o,i_o,RAAN_o,u_o, ...
    d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y)
    
    [a_t,e_x_t,e_y_t,i_t,RAAN_t,u_t] = ...
    ROE2quasi_nonsing(a_init,ex_init,ey_init,i_init,RAAN_init,u_init, ...
    d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y);
    
    [a_t, e_t, i_t, RAAN_t, w_t, nu_t, M_t] = ...
    quasi_nonsing2OE(a_t, e_x_t, e_y_t, i_t, RAAN_t, u_t);
    
    [r_ECI, v_ECI] = OE2ECI(a_t, e_t, i_t, RAAN_t, w_t, nu_t);
end