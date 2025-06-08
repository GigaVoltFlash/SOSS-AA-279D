function [x_update, P_update, y_pred, y_post] = ekf_roes_with_control(x_prior, y_meas, P_prior, SV1_state, SV1_OE_state, Q, R, dt, u)

    a_o = SV1_OE_state(1);
    e_x_o = SV1_OE_state(2);
    e_y_o = SV1_OE_state(3);
    i_o = SV1_OE_state(4);
    RAAN_o = SV1_OE_state(5);
    u_o = SV1_OE_state(6);

    B = Steindorf_control_input_matrix(SV1_OE_state);

    [a_o,e_o,i_o,RAAN_o,w_o,nu_o, ~] = quasi_nonsing2OE(a_o, e_x_o, e_y_o, i_o, RAAN_o, u_o); 
    
    SV1_OE_sing = [a_o,e_o,i_o,RAAN_o,w_o,nu_o];

    [STM_big,~] = roe_stm_j2(dt, x_prior', SV1_OE_sing);
    STM_curr = squeeze(STM_big);
    x_bar = STM_curr*x_prior + B*u';

    % EKF Covariance Prediction
    P_bar = STM_curr * P_prior * STM_curr' + Q;

    % EKF Predict Measurements
    y_pred = measurement_model_SV3(x_bar,SV1_OE_state);

    % EKF Update
    H = measurement_sensitivity_matrix_SV3(SV1_state(1:3),SV1_state(4:6),SV1_OE_state);
    K = (P_bar*(H'))/(H*P_bar*H' + R);

    x_update = x_bar+(K*(y_meas-y_pred));
    P_update = (eye(6)-K*H)*P_bar*(eye(6)-K*H)' + K*R*K';
    y_post = measurement_model_SV3(x_update,SV1_OE_state);
end