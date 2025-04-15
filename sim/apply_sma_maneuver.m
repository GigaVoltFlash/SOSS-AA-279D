function [rho_RTN_w_maneuver, rho_RTN_dot_w_maneuver, t_combined] = apply_sma_maneuver(state_abs_init,tstart,tint,maneuver_time,tend, ...
    desired_delta_a,r_ECI_no_j2,v_ECI_no_j2)
    [t_1, state1] = rk4_eom_ECI(tstart:tint:(maneuver_time-tint), state_abs_init, false); % no J2
    r_ECI_1 = state1(:,1:3);
    v_ECI_1 = state1(:,4:6);
    
    [a_1,e_1,i_1,RAAN_1,omega_1,nu_1] = ECI2OE_array(r_ECI_1,v_ECI_1);
    periapsis_indices = find(nu_1 < 0.5); 
    last_periapsis_idx = periapsis_indices(end);
    
    r_inter = state1(last_periapsis_idx,1:3);
    v_inter = state1(last_periapsis_idx,4:6);
    t_inter = t_1(last_periapsis_idx);
    t_first = t_1(1:last_periapsis_idx);
    state_first = state1(1:last_periapsis_idx,:);
    
    [a_inter,e_inter] = ECI2OE(r_inter,v_inter);
    delta_v_t = sma_maneuver(a_inter,e_inter,r_inter,desired_delta_a); 
    delta_v_RTN = [0;delta_v_t;0];
    delta_v_ECI = dv_RTN2ECI(r_inter, v_inter, delta_v_RTN);
    
    v_inter_new = v_inter+delta_v_ECI';
    state_abs_inter_new = [r_inter';v_inter_new'];
    
    [t_2, state2] = rk4_eom_ECI(t_inter+tint:tint:tend, state_abs_inter_new, false);
    
    state_combined = [state_first; state2];  
    t_combined = [t_first,t_2]; 
    
    r_ECI_new = state_combined(:,1:3);
    v_ECI_new = state_combined(:,4:6);
    
    [rho_RTN_w_maneuver, rho_RTN_dot_w_maneuver] = ECI2RTN_rel(r_ECI_no_j2, v_ECI_no_j2, r_ECI_new, v_ECI_new);
end