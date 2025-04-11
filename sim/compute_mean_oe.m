function [a_3,e_3,i_3,RAAN_3,omega_3,nu_3,t_3] = compute_mean_oe(a_init,ex_init,ey_init,i_init,RAAN_init,u_init,tstart,tint,tend)
    quasi_nonsing_oe_init = [a_init;ex_init;ey_init;i_init;RAAN_init;u_init];
    
    [t_3, state3] = rk4_secular_J2_GVE(tstart:tint:tend, quasi_nonsing_oe_init);
    
    a_3 = state3(:,1);
    e_x_3 = state3(:,2);
    e_y_3 = state3(:,3);
    i_3 = state3(:,4);
    RAAN_3 = state3(:,5);
    u_3 = state3(:,6);
    
    e_3 = sqrt(e_x_3.^2 + e_y_3.^2);             
    omega_3 = rad2deg(atan2(e_y_3, e_x_3));  
    
    nu_3 = u_3 - omega_3;
    nu_3 = mod(nu_3, 360);
