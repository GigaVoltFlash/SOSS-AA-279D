%% AA279D Control for Mode 1: Station Keeping
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

% Mode 3 control is for the Docker (SV3)'s final docking procedure towards SV1
% Since this mode is a small maneuver, we use a closed-form solution with radial delta v inputs to get the 
% the desired control. The closed-form solution gives us the optimal times to perform this maneuver as well.
function [delta_v_vals, delta_v_times] = mode3_control(roe_initial, roe_final, SV1_OE_init)

    % USE CLOSED-FORM SOLUTION
    t0 = 0;
    tf = 100; % No idea, just putting something for now.
    [delta_v_vals, delta_v_times] = closed_form_control(roe_initial, roe_final, SV1_OE_init, t0, tf);


end