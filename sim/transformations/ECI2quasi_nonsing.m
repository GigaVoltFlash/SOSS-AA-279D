function OE_qns = ECI2quasi_nonsing(ECI_state)
     r = ECI_state(1:3);
     v = ECI_state(4:6);

    [a, e, i, RAAN, omega, nu, M] = ECI2OE(r,v);

    [a, e_x, e_y, i, RAAN, u] = OE2quasi_nonsing(a, e, i, RAAN, omega, M);

    OE_qns = [a, e_x, e_y, i, RAAN, u]; 
end