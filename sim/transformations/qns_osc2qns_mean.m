function mean_qns_elem = qns_osc2qns_mean(osc_qns_elem)
    % Converts osculating qns OE to mean qns OE
    a_osc = osc_qns_elem(1);
    ex_osc = osc_qns_elem(2);
    ey_osc = osc_qns_elem(3);
    i_osc = osc_qns_elem(4);
    RAAN_osc = osc_qns_elem(5);
    u_osc = osc_qns_elem(6);

    [a_osc,e,i_osc,RAAN_osc,w_osc,nu_osc, M_deg] = quasi_nonsing2OE(a_osc, ex_osc, ey_osc, i_osc, RAAN_osc, u_osc);

    osc_elem = [a_osc*1e3,e,deg2rad(i_osc),deg2rad(RAAN_osc),deg2rad(w_osc),deg2rad(M_deg)];

    J2_flag = true;
    mean_elem = osc2mean(osc_elem,J2_flag);

    a_mean = mean_elem(1)/1e3; % m --> km
    e_mean = mean_elem(2);
    i_mean = mean_elem(3); % rad
    RAAN_mean = mean_elem(4); % rad
    w_mean = mean_elem(5); % rad
    M_mean = mean_elem(6); % rad

    [a_mean, ex_mean, ey_mean, i_mean, RAAN_mean, u_mean] = OE2quasi_nonsing(a_mean, e_mean, rad2deg(i_mean), rad2deg(RAAN_mean), rad2deg(w_mean), rad2deg(M_mean));

    mean_qns_elem = [a_mean, ex_mean, ey_mean, i_mean, RAAN_mean, u_mean];
end