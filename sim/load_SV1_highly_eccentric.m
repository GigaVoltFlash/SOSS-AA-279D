%% SV1 Initial Highly Eccentric Conditions
% Loads the initial conditions of the chief satellite, i.e. SV1
% Units: kilometers and degrees

a_SV1_init = 6944;
ex_SV1_init = -4e-5;
ey_SV1_init = 0.55; %0.15;
i_SV1_init = 99.4;
RAAN_SV1_init = -151.1;
u_SV1_init = -47.9;
[a_SV1_init, e_SV1_init, i_SV1_init, RAAN_SV1_init, w_SV1_init, nu_SV1_init, M_SV1_init] = quasi_nonsing2OE(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init);

[r_SV1_ECI_init, v_SV1_ECI_init] = OE2ECI(a_SV1_init, e_SV1_init, i_SV1_init, RAAN_SV1_init, w_SV1_init, nu_SV1_init);
state_init = [r_SV1_ECI_init; v_SV1_ECI_init];
