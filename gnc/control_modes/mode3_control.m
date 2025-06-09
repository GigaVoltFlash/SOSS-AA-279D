%% AA279D Control for Mode 3: Going to proximity
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

% Mode 3 control is for the Docker (SV3)'s final docking procedure towards SV1
% Since this mode is a small maneuver, we use a closed-form solution with tangential v inputs to get the 
% the desired control. The closed-form solution gives us the optimal times to perform this maneuver as well.
function [delta_v_vals, delta_v_times] = mode3_control(roe_initial, roe_final, init_time, final_time, SV1_OE_init, u_SV1_init)
    % n_man = 3; % number of maneuvers we want to do
    % 
    % % Picking the control actions relative to the initial time
    % delta_v_times = linspace(init_time, final_time, n_man+1); % Should give me equally spaced time-steps from start to end times?
    % delta_v_times = delta_v_times(2:end);
    % delta_v_vals = naive_least_squares(delta_v_times, roe_initial, roe_final, SV1_OE_init, u_SV1_init, init_time, final_time);

    % Another method: Use D'Amico's formulation for tangential and normal burns. This one is more optimal but takes longer
    STM_from_init = calc_STM_for_control(init_time, final_time, SV1_OE_init);
    no_control_roe_final = STM_from_init * roe_initial';
    total_delta_roe = roe_final' - no_control_roe_final; % Total roe change
    [delta_v_vals_inplane, delta_v_times_inplane] = triplet_burn_tangential(total_delta_roe(1), total_delta_roe(3:4), total_delta_roe(2), SV1_OE_init, u_SV1_init, init_time, final_time);
    [delta_v_vals_outplane, delta_v_times_outplane] = single_impulse_crosstrack(total_delta_roe(5:6), SV1_OE_init, u_SV1_init, init_time, final_time);
    delta_v_vals = [delta_v_vals_inplane; delta_v_vals_outplane];
    delta_v_times = [delta_v_times_inplane, delta_v_times_outplane];
end
