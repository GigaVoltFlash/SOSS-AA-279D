function delta_v_vals_sk = station_keeping_continuous(SV_roe, SV_roe_nom, SV_delta_de_max, SV_delta_di_max, SV1_oe)
    % Calculate offsets for station keeping SV
    SV_dphi = asin(SV_delta_de_max/norm(SV_roe_nom(3:4)));
    SV_de_des = [SV_roe_nom(3)*cos(SV_dphi) - SV_roe_nom(4)*sin(SV_dphi);...
            SV_roe_nom(3)*sin(SV_dphi) + SV_roe_nom(4)*cos(SV_dphi)];

    delta_v_vals_ecc = zeros(3,1);
    delta_v_vals_inc = zeros(3,1);
    % SV eccentricity station keeping
    if norm(SV_roe(3:4) - SV_roe_nom(3:4)) > SV_delta_de_max
        SV_des_roe = SV_roe_nom;
        SV_des_roe(3:4) = SV_de_des;
        % perform station keeping
        delta_v_vals_ecc = Lyapunov_feedback_control(SV_roe, SV_des_roe, SV1_oe, 4, 1000);
    end

    % SV inclination station keeping
    if norm(SV_roe(5:6) - SV_roe_nom(5:6)) > SV_delta_di_max
        SV_des_roe = SV_roe_nom;
        SV_di_des = [SV_roe_nom(5);SV_roe_nom(6) - sign(SV_roe(5))*SV_delta_di_max];
        SV_des_roe(5:6) = SV_di_des;
        % perform station keeping 
        delta_v_vals_inc = Lyapunov_feedback_control(SV_roe, SV_des_roe, SV1_oe, 4, 1000);
    end

    delta_v_vals_sk = delta_v_vals_ecc + delta_v_vals_inc; % Make sure that this is ok, and that ecc control is tangential and inc control is normal
end