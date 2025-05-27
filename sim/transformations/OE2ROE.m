function roe = OE2ROE(chief_oe_qns, deputy_oe_qns)

    a_o = chief_oe_qns(1);
    e_x_o = chief_oe_qns(2);
    e_y_o = chief_oe_qns(3);
    i_o = chief_oe_qns(4);
    RAAN_o = chief_oe_qns(5);
    u_o = chief_oe_qns(6);
    a_t = deputy_oe_qns(1);
    e_x_t = deputy_oe_qns(2);
    e_y_t = deputy_oe_qns(3);
    i_t = deputy_oe_qns(4);
    RAAN_t = deputy_oe_qns(5);
    u_t = deputy_oe_qns(6);


    [d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y] = quasi_nonsing2ROE(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o, ...
    a_t, e_x_t, e_y_t, i_t, RAAN_t, u_t);

    roe = [d_a,d_lambda,d_e_x,d_e_y,d_i_x,d_i_y];
end