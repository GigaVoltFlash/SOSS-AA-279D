%% AA279D Control for Mode 2: Approach
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

% Mode 2 control is for the Docker (SV3)'s approach towards SV1
% Since this mode is a large maneuver, we use naive least squares to break it down into smaller maneuvers
% Each of these maneuvers is then achieved using s
function [delta_v_vals, delta_v_times] = mode2_control(roe_initial, roe_final, init_time, final_time, SV1_OE_init)
    % Need to see whether to take initial time or initial mean argument of latitude.
    % The output will specifi

    % This is going to be an along-track burn to get the desired effect (no idea what this is actually going to do)
    % This should give specific time-steps 

    % Steps here:
    % 1. ####### Calculate the guidance law to get from start to finish #######
    n_man = 3; % number of maneuvers we want to do

    % Picking the guidance law times relative to the initial time
    delta_v_times = linspace(init_time, final_time, n_man+1); % Should give me equally spaced time-steps from start to end times?
    delta_v_times = delta_v_times(2:end);
    delta_v_vals = naive_least_squares(delta_v_times, roe_initial, roe_final, SV1_OE_init, init_time, final_time);

end

