% Tried a bunch of things because I was getting only a small radial
% component in the RTN velocity. But I realized that since the RTN frame is
% rotating with the satellite, there is no tangential velocity in that
% direction.

% function [r_RTN, v_RTN]= ECI2RTN(r_ECI,v_ECI)
%     % Takes in r_ECI and v_ECI has time arrays of shape nx3
%     n = size(r_ECI, 1);
%     R_hat = r_ECI./vecnorm(r_ECI, 2, 2);
%     N = cross(r_ECI, v_ECI, 2);
%     N_hat = N./vecnorm(N, 2, 2);
%     T_hat = cross(N_hat, R_hat, 2);
% 
%     % SatelliteDynamics.jl from SISL
%     fdot = vecnorm(cross(r_ECI, v_ECI, 2), 2, 2) ./ (vecnorm(r_ECI, 2, 2)).^2;
%     omega = horzcat(zeros(n, 1), zeros(n, 1), fdot);
% 
%     Q_eci2rtn = zeros(n, 3, 3); 
%     Q_eci2rtn(:,1,:) = R_hat; 
%     Q_eci2rtn(:,2,:) = T_hat; 
%     Q_eci2rtn(:,3,:) = N_hat;
%     % Q_eci2rtn(:, :, 1) = R_hat;  % Column 1
%     % Q_eci2rtn(:, :, 2) = T_hat;  % Column 2
%     % Q_eci2rtn(:, :, 3) = N_hat;  % Column 3
% 
%     r_RTN = zeros(n, 3);
%     v_RTN = zeros(n, 3);
%     for i = 1:n
%         r_RTN(i,:) = squeeze(Q_eci2rtn(i, :, :)) * r_ECI(i, :)';
%         % v_RTN(i,:) = squeeze(Q_eci2rtn(i, :, :)) * (v_ECI(i, :)' - cross(omega(i, :), r_ECI(i,:))');
%         v_RTN(i,:) = squeeze(Q_eci2rtn(i, :, :)) * v_ECI(i, :)' - cross(squeeze(Q_eci2rtn(i, :, :))' * omega(i, :)', r_RTN(i,:))'; % mathematically correct one
%     end
% end
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
    omega_eci = cross(r_ECI, v_ECI, 2) ./ vecnorm(r_ECI, 2, 2).^2;  % each row is ω_eci (in ECI)

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

