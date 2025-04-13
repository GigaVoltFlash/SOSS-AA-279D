
function [rho_RTN, rho_RTN_dot] = ECI2RTN_rel_init(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1)
    % _0 indicates input position and velocities of the chief
    % _1 indicates input position and velocities of the deputy

    % The function returns the relative position and velocity of the deputy
    % in the chief's 

    % Build RTN basis vectors
    R_hat = r_ECI_0 ./ norm(r_ECI_0);
    h_vec = cross(r_ECI_0, v_ECI_0);
    N_hat = h_vec ./ norm(h_vec);
    T_hat = cross(N_hat, R_hat);

    Q_eci2rtn = [R_hat', T_hat', N_hat']';

    % Angular velocity of RTN frame (in ECI frame)
    omega0_eci = cross(r_ECI_0, v_ECI_0) ./ norm(r_ECI_0).^2; 

    rho_ECI = r_ECI_1 - r_ECI_0; 

    rho_RTN = Q_eci2rtn * rho_ECI';
    rho_dot_corrected = (v_ECI_1 - v_ECI_0)' - cross(omega0_eci', r_ECI_0');
    rho_RTN_dot = Q_eci2rtn * rho_dot_corrected;
end