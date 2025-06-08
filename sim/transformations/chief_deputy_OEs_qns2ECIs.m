function [SV1_state, SV2_state, SV3_state, ...
          r_ECI_SV1, v_ECI_SV1, r_ECI_SV2, v_ECI_SV2, r_ECI_SV3, v_ECI_SV3] = ...
          chief_deputy_OEs_qns2ECIs(SV1_OE_state, SV2_OE_state, SV3_OE_state)
% CHIEF_DEPUTY_OES_QNS2ECIS Converts chief and deputy QNS orbital elements to ECI states
%
% Inputs:
%   SV1_OE_state - 1x6 vector [a [km], ex, ey, i [deg], RAAN [deg], u [deg]] chief
%   SV2_OE_state - 1x6 vector [a [km], ex, ey, i [deg], RAAN [deg], u [deg]] deputy 2
%   SV3_OE_state - 1x6 vector [a [km], ex, ey, i [deg], RAAN [deg], u [deg]] deputy 3
%
% Outputs:
%   SV1_state - 6x1 vector [r_ECI; v_ECI] of chief (position and velocity in km and km/s)
%   SV2_state - 6x1 vector [r_ECI; v_ECI] of deputy 2
%   SV3_state - 6x1 vector [r_ECI; v_ECI] of deputy 3
%   r_ECI_SV1, v_ECI_SV1 - 3x1 vectors of chief position and velocity
%   r_ECI_SV2, v_ECI_SV2 - 3x1 vectors of deputy 2 position and velocity
%   r_ECI_SV3, v_ECI_SV3 - 3x1 vectors of deputy 3 position and velocity

    % --- Chief ---
    a_o    = SV1_OE_state(1);
    e_x_o  = SV1_OE_state(2);
    e_y_o  = SV1_OE_state(3);
    i_o    = SV1_OE_state(4);
    RAAN_o = SV1_OE_state(5);
    u_o    = SV1_OE_state(6);

    % --- Deputy 2 ---
    a_2    = SV2_OE_state(1);
    e_x_2  = SV2_OE_state(2);
    e_y_2  = SV2_OE_state(3);
    i_2    = SV2_OE_state(4);
    RAAN_2 = SV2_OE_state(5);
    u_2    = SV2_OE_state(6);

    % --- Deputy 3 ---
    a_3    = SV3_OE_state(1);
    e_x_3  = SV3_OE_state(2);
    e_y_3  = SV3_OE_state(3);
    i_3    = SV3_OE_state(4);
    RAAN_3 = SV3_OE_state(5);
    u_3    = SV3_OE_state(6);

    % --- Chief ---
    [a_o,e_o,i_o,RAAN_o,w_o,nu_o, ~] = quasi_nonsing2OE(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o); 
    [r_ECI_SV1,v_ECI_SV1] = OE2ECI(a_o,e_o,i_o,RAAN_o,w_o,nu_o);
    
    % --- Deputy 2 ---
    [a_2,e_2,i_2,RAAN_2,w_2,nu_2, ~] = quasi_nonsing2OE(a_2, e_x_2, e_y_2, i_2, RAAN_2, u_2); 
    [r_ECI_SV2,v_ECI_SV2] = OE2ECI(a_2,e_2,i_2,RAAN_2,w_2,nu_2);

    % --- Deputy 3 ---
    [a_3,e_3,i_3,RAAN_3,w_3,nu_3, ~] = quasi_nonsing2OE(a_3, e_x_3, e_y_3, i_3, RAAN_3, u_3); 
    [r_ECI_SV3,v_ECI_SV3] = OE2ECI(a_3,e_3,i_3,RAAN_3,w_3,nu_3);

    % --- Output states ---
    SV1_state = [r_ECI_SV1' , v_ECI_SV1'];
    SV2_state = [r_ECI_SV2' , v_ECI_SV2'];
    SV3_state = [r_ECI_SV3' , v_ECI_SV3'];

end
