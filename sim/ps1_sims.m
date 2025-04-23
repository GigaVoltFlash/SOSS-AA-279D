%% AA279D Simulations for PSET 1
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

%%%%%%%%%%% RUN SIM WITH AND WITHOUT J2 %%%%%%%%%%%
[t_1, state1] = rk4_eom_ECI(tstart:tint:tend, state_init, false);
[t_2, state2] = rk4_eom_ECI(tstart:tint:tend, state_init, true);
r_ECI_with_j2 = state2(:,1:3);
v_ECI_with_j2 = state2(:,4:6);
r_ECI_no_j2 = state1(:,1:3);
v_ECI_no_j2 = state1(:,4:6);

%%%%%%% KEPLERIAN PROPAGATION %%%%%%% 
% Propagates over the same times t_1 as the sim without J2 
[r_ECI_keplerian, v_ECI_keplerian] = keplerian_propagator(a_SV1_init, e_SV1_init, i_SV1_init, RAAN_SV1_init, w_SV1_init, t_1, M_SV1_init);

%%%%%% CONVERT TO RTN %%%%%%%%
[r_RTN_no_j2, v_RTN_no_j2] = ECI2RTN(r_ECI_no_j2,v_ECI_no_j2);
[r_RTN_keplerian, v_RTN_keplerian] = ECI2RTN(r_ECI_keplerian,v_ECI_keplerian);