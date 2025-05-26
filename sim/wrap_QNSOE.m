function QNSOE_wrapped = wrap_QNSOE(QNSOE_deg)
% wrap_QNSOE Wraps the u (argument of latitude) in quasi-nonsingular OE to [0, 360)
% Input:
%   QNSOE_deg - Nx6 matrix or 1x6 vector: [a, e_x, e_y, i_x, i_y, u] in degrees
% Output:
%   QNSOE_wrapped - same shape as input, with u wrapped to [0, 360)

    % Determine if input is a row vector or matrix
    is_row_vector = isvector(QNSOE_deg) && size(QNSOE_deg,1) == 1;

    % Copy input
    QNSOE_wrapped = QNSOE_deg;

    % Wrap only the 6th element (u)
    QNSOE_wrapped(:,6) = mod(QNSOE_deg(:,6), 360);

    % Preserve row vector shape if applicable
    if is_row_vector
        QNSOE_wrapped = QNSOE_wrapped(1,:);
    end
end