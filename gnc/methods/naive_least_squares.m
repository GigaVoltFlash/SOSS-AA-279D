
% For a given set of times, this function uses naive least squares to find a set of delta v maneuvers 
% that brings you close to the change in ROE that is desired. We can try and use this for the approach, or for
% a reconfiguration
function [delta_v_vals] = naive_least_squares(t_maneuvers, roe_init, roe_final, SV1_oe_init, t_final)
    m = length(t_maneuvers);
    total_delta_roe = (roe_final - roe_init)'; % Total roe change
    Matrix_block = zeros(6, 3*m);

    for i=1:m
        t_val = t_maneuvers(i);
        STM = calc_STM_for_control(t_val, t_final, SV1_oe_init);
        Gamma = calc_Gamma_for_control(t_val, SV1_oe_init);
        Matrix_block(:, 3*(i-1)+1:3*i) = STM * Gamma;
    end
 
    delta_v_vals = reshape(pinv(Matrix_block) * total_delta_roe,m,3);
    
end

% Old version that did waypoints. We can also do this, and use closed form to find the in-between maneuvers.
% function [roe_waypoints] = naive_least_squares(t_maneuvers, roe_init, roe_final, SV1_oe_init)
   
%     m = length(t_maneuvers);
%     % SV1 oe vals is another 5 vector of the orbital elements of the
%     % chief at the initial time
%     total_delta_roe = roe_final - roe_init; % Total roe change
%     % to solve for 6m variables, we use least squares
%     STM_block = zeros(6, 6*m);
%     [STMs, expected_roe_states_garbage] = roe_stm_j2(t_maneuvers, roe_init, SV1_oe_init);
    
%     % % Re-structuring the STM to be 6x6m as specified in the slides
%     % for i=1:m
%     %     for j=1:6
%     %         STM_block(:, m*(j-1) + i) = STMs(i, :, j);
%     %     end
%     % end
    
%     % More efficient re-structuring from ChatGPT
%     STM_block = permute(STMs, [3, 2, 1]);
%     STM_block = reshape(STM_block, 6, []);

%     % Least squares analytical solution
%     roe_waypoints = pinv(STM_block) * total_delta_roe;
% end