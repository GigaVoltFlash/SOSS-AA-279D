%% AA279D Control for Mode 2: Approach
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

function [delta_v_vals, delta_v_times] = mode2_control(roe_initial, roe_final, init_time, final_time, SV1_OE_init, u_SV1_init)
    n_man = 3; % number of maneuvers we want to do

    % Picking the control actions relative to the initial time
    delta_v_times = linspace(init_time, final_time, n_man+1); % Should give me equally spaced time-steps from start to end times?
    delta_v_times = delta_v_times(2:end);
    delta_v_vals = naive_least_squares(delta_v_times, roe_initial, roe_final, SV1_OE_init, u_SV1_init, init_time, final_time);

    % % Another method: Use D'Amico's formulation for tangential and normal burns. This one is more optimal but takes longer
    % roe_diff = roe_final - roe_initial;
    % [delta_v_vals, delta_v_times] = damico_maneuvers(roe_diff(1), roe_diff(3:4), roe_diff(5:6), SV1_OE_init, u_SV1_init, init_time, final_time);
end

