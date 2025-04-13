function [rho_RTN, rho_RTN_dot] = ECI2RTN_rel(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1)
    % _0 indicates input position and velocities of the chief
    % _1 indicates input position and velocities of the deputy

    % The function returns the relative position and velocity of the deputy
    % in the chief's 
    n = size(r_ECI_0, 1);

    % Build RTN basis vectors
    R_hat = r_ECI_0 ./ vecnorm(r_ECI_0, 2, 2);
    h_vec = cross(r_ECI_0, v_ECI_0, 2);
    N_hat = h_vec ./ vecnorm(h_vec, 2, 2);
    T_hat = cross(N_hat, R_hat, 2);

    Q_eci2rtn = zeros(n, 3, 3);
    for i = 1:n
        Q_eci2rtn(i, :, :) = [R_hat(i,:)', T_hat(i,:)', N_hat(i,:)']';
    end

    % Angular velocity of RTN frame (in ECI frame)
    omega0_eci = cross(r_ECI_0, v_ECI_0, 2) ./ vecnorm(r_ECI_0, 2, 2).^2;  % each row is ω_eci (in ECI)

    rho_ECI = r_ECI_1 - r_ECI_0; 

    % Output
    rho_RTN = zeros(n, 3);
    rho_RTN_dot = zeros(n, 3);

    for i = 1:n
        Q = squeeze(Q_eci2rtn(i, :, :));   % ECI → RTN rotation
        rho_RTN(i, :) = Q * rho_ECI(i, :)';
        rho_dot_corrected = (v_ECI_1(i,:) - v_ECI_0(i,:))' - cross(omega0_eci(i,:)', r_ECI_0');
        rho_RTN_dot(i,:) = Q * rho_dot_corrected;
    end
end
