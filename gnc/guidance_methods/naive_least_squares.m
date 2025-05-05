function [roe_waypoints] = naive_least_squares(t_maneuvers, roe_init, roe_final, SV1_oe_init)
    
   
    m = length(t_maneuvers);
    % SV1 oe vals is another 5 vector of the orbital elements of the
    % chief at the initial time
    total_delta_roe = roe_final - roe_init; % Total roe change
    % to solve for 6m variables, we use least squares
    STM_block = zeros(6, 6*m);
    [STMs, expected_roe_states_garbage] = roe_stm_j2(t_maneuvers, roe_init, SV1_oe_init);
    
    % % Re-structuring the STM to be 6x6m as specified in the slides
    % for i=1:m
    %     for j=1:6
    %         STM_block(:, m*(j-1) + i) = STMs(i, :, j);
    %     end
    % end
    
    % More efficient re-structuring from ChatGPT
    STM_block = permute(STMs, [3, 2, 1]);
    STM_block = reshape(STM_block, 6, []);

    % Least squares analytical solution
    roe_waypoints = pinv(STM_block) * total_delta_roe;
end