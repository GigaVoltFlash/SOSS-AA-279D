function [r_RTN, v_RTN] = ECI2RTN(r_ECI, v_ECI)
    n = size(r_ECI, 1);

    % Build RTN basis vectors
    R_hat = r_ECI ./ vecnorm(r_ECI, 2, 2);
    h_vec = cross(r_ECI, v_ECI, 2);
    N_hat = h_vec ./ vecnorm(h_vec, 2, 2);
    T_hat = cross(N_hat, R_hat, 2);

    Q_eci2rtn = zeros(n, 3, 3);
    for i = 1:n
        Q_eci2rtn(i, :, :) = [R_hat(i,:)', T_hat(i,:)', N_hat(i,:)']';
    end

    % Angular velocity of RTN frame (in ECI frame)
    omega_eci = cross(r_ECI, v_ECI, 2) ./ vecnorm(r_ECI, 2, 2).^2;

    % Output
    r_RTN = zeros(n, 3);
    v_RTN = zeros(n, 3);

    for i = 1:n
        Q = squeeze(Q_eci2rtn(i, :, :));   % ECI → RTN rotation
        r_RTN(i,:) = Q * r_ECI(i, :)';
        v_corrected = v_ECI(i,:)' - cross(omega_eci(i,:)', r_ECI(i,:)');
        v_RTN(i,:) = Q * v_corrected;
    end
end

