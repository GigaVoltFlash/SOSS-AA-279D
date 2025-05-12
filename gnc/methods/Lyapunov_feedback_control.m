function control_vec = Lyapunov_feedback_control(A_matrix, ROE, BINV, ERROR, OE_CHIEF, PFEED, myDrag)
    % Inputs:
    % A_matrix: 5x5 matrix
    % ROE: 5x1 vector
    % BINV: 2x5 matrix (stored column-wise)
    % ERROR: 5x1 vector
    % OE_CHIEF: 6x1 vector
    % PFEED: 5x5 matrix
    % myDrag: 5x1 vector
    
    % Output:
    % out_dep: 3x1 control output

    % Extract scalar normalization factor
    oe_norm = OE_CHIEF(1);  % Corresponds to oe_chief(0) in C

    % Normalize inputs
    roe = ROE ./ oe_norm;
    drag = myDrag ./ oe_norm;
    error = ERROR ./ oe_norm;

    % Extract matrices
    Aplant = reshape(A_matrix, 5, 5);     % Assuming row-major flattening
    P = reshape(PFEED, 5, 5);
    Binv = reshape(BINV, 2, 5);

    % Compute control input
    u_2 = -Binv * (Aplant * roe + drag + P * error);

    % Output
    control_vec = [0; u_2(:)];  % 3x1 vector: [0; u_2(1); u_2(2)]
end
