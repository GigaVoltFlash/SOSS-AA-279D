%% AA279D
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

% This file has the implementation of a closed-form control solution, that is used extensively in the proximity
% operations and potentially in the approach calculations.
% ANSHUK CODED VERSION. SEE PROVIDED METHOD INSTEAD.

function [delta_v_vals, delta_v_times] = closed_form_control(roe_initial, roe_final, SV1_oe_init, t0, tf)
    % Some relevant parameters that are useful
    global J2 R_earth mu_earth
    a = SV1_oe_init(1);
    ecc = SV1_oe_init(2);
    i_SV1_init = SV1_oe_init(3);
    i = deg2rad(i_SV1_init);
    w_SV1_init = SV1_oe_init(5);
    w = deg2rad(w_SV1_init);
    eta = sqrt(1 - ecc^2);
    kappa = 3/4*J2*R_earth^2 * sqrt(mu_earth)/(a^(7/2) * eta^4);
    Q = 5*(cos(i))^2 - 1;
    n = sqrt(mu_earth/a^3);
    omega_dot = kappa * Q;
    
    % Step 1. Identify dominant case 
    STM_from_init = calc_STM_for_control(t0, t_final, SV1_oe_init);
    no_control_roe_final = STM_from_init * roe_init';
    roe_diff = roe_final - no_control_roe_final;
    d_delta_a = roe_diff(1);
    d_delta_lambda = roe_diff(2);
    d_delta_e = roe_diff(3:4);
    d_delta_i = roe_diff(5:6);

    a_dv_min = abs(a*d_delta_a)*n/2; % From Table 5.13 in Chernick
    m = -2*abs(d_delta_a0)/abs(d_delta_lambda0);
    lambda_dv_min = n*a*abs((m*d_delta_lambda - d_delta_a)/(d_delta_a0)); % What's delta M here?
    ecc_dv_min = norm(a*d_delta_e)*n/2;
    i_dv_min = norm(a*d_delta_i)*n; % Do this separately but calculate it in case we need a cross-track burn

    dv_lb = max(a_dv_min, lambda_dv_min, ecc_dv_min); % Get the max and the argmax
    dv_lb_case = argmax(a_dv_min, lambda_dv_min, ecc_dv_min);

    % Based on the dominant case, approach the rest
    if dv_lb_case==1
        % Semi-major axis is dominant
        m = 3; % Number of maneuvers required
        if d_delta_a > 0
            sign = 1;
        else
            sign = -1;
        end
        dv_vec = [0, sign*a_dv_min, 0]; % Tangential is optimal
        dv_vec = repmat(dv_vec, m, 1);
        offset = ceiling((-atan(d_delta_e(2)/d_delta_e(1)) + (omega_dot*tf + w))/pi); % Calculate a pi that works
        T_opt_vals = [];     % Initialize list of valid T values
        k_vals = [];    % Corresponding k values
        k = 0;
        while true
            T_k = ((atan(d_delta_e(2)/d_delta_e(1)) + offset*pi) - (omega_dot*tf + w))/omega_dot + k*pi/omega_dot;
            if T_k > tf
                break;
            end
            T_opt_vals(end+1) = T_k;
            k_vals(end+1) = k;
            k = k + 1;
        end
        num_options = length(T_opt_vals);
        dv_vec = repmat(dv_vec, num_options, 1);
        sets = zeros(num_options, 6);
        % Create reachable set
        for iter=1:num_options
            T_opt = T_opt_vals(iter);
            STM = calc_STM_for_control(T_opt, tf, SV1_oe_init);
            Gamma = calc_Gamma_for_control(T_opt, SV1_oe_init);
            sets(iter, :) = STM * Gamma * dv_vec(iter, :);
        end

        % Try all combinations 
        combo_idx = nchoosek(1:num_options, m);   % all n choose m combinations of row indices
        num_combos = size(combo_idx, 1);
        good_c = [0, 0, 0];
        good_combo = [0, 0, 0];
        for i = 1:num_combos
            selected_rows = combo_idx(i, :);   % row indices for this combination
            subset = sets(selected_rows, :);   % m x 6 matrix of selected vectors
            A = [1, 1, 1;...
                 a*subset(1, 2), a*subset(2, 2), a*subset(3, 2);...
                 sign*a*norm(subset(1, 3:4)), sign*a*norm(subset(2, 3:4)), sign*a*norm(subset(3, 3:4))];

            b = [1; a*d_delta_lambda; a*norm(d_delta_e)];

            c = A\b; % Get the maneuver magnitudes

            if all(c) > 0 && sum(c) == 1
                good_c = c;
                good_combo = selected_rows;
                break
            end
        end

        delta_v_vals = dv_vec(good_combo, :) .* good_c(:);% TODO: Make sure this is multiplying correctly
        delta_v_times = T_opt_vals(good_combo);

    elseif d_lb_case==2
        % ADD DELTA DELTA LAMBDA CASE

    else
        % Delta ecc is dominant
        m = 3; % Number of maneuvers required
        offset = ceiling((-atan(d_delta_e(2)/d_delta_e(1)) + (omega_dot*tf + w))/pi); % Calculate a pi that works
        T_opt_vals = [];     % Initialize list of valid T values
        k_vals = [];    % Corresponding k values
        k = 0;
        while true
            T_k = ((atan(d_delta_e(2)/d_delta_e(1)) + offset*pi) - (omega_dot*tf + w))/omega_dot + k*pi/omega_dot;
            if T_k > tf
                break;
            end
            T_opt_vals(end+1) = T_k;
            k_vals(end+1) = k;
            k = k + 1;
        end
        num_options = length(T_opt_vals);
        dv_vec = zeros(num_options, 3);
        for j=1:num_options
            if rem(k_vals(j) + offset, 2) == 0
                dv_vec(j, :) = [0, ecc_dv_min, 0]; % Tangential is optimal
            else
                dv_vec(j, :) = [0, -ecc_dv_min, 0]; % Tangential is optimal
            end
        end
        
        sets = zeros(num_options, 6);
        % Create reachable set
        for iter=1:num_options
            T_opt = T_opt_vals(iter);
            STM = calc_STM_for_control(T_opt, tf, SV1_oe_init);
            Gamma = calc_Gamma_for_control(T_opt, SV1_oe_init);
            sets(iter, :) = STM * Gamma * dv_vec(iter, :);
        end

        % Try all combinations 
        combo_idx = nchoosek(1:num_options, m);   % all n choose m combinations of row indices
        num_combos = size(combo_idx, 1);
        good_c = [0, 0, 0];
        good_combo = [0, 0, 0];
        for i = 1:num_combos
            selected_rows = combo_idx(i, :);   % row indices for this combination
            subset = sets(selected_rows, :);   % m x 6 matrix of selected vectors
            A = [1, 1, 1;...
                 a*subset(1, 2), a*subset(2, 2), a*subset(3, 2);...
                 a*subset(1, 1), a*subset(2, 1), a*subset(3, 1)];

            b = [1; a*d_delta_lambda; a*d_delta_a];

            c = A\b; % Get the maneuver magnitudes

            if all(c) > 0 && sum(c) == 1
                good_c = c;
                good_combo = selected_rows;
                break
            end
        end

        delta_v_vals = dv_vec(good_combo, :) .* good_c(:);% TODO: Make sure this is multiplying correctly
        delta_v_times = T_opt_vals(good_combo);
    end

    % Find the optimal time for inclination change maneuver
    % TODO: ADD INCLINATION CHANGE STUFF HERE
    if abs(d_delta_i) > 1e-4 % Or whatever relative err we want to deal with here to do a cross-track maneuver
        offset = ceiling((-atan(d_delta_i(2)/d_delta_i(1)) + w)/pi); % Calculate a pi that works
        T_opt_vals = [];     % Initialize list of valid T values
        k_vals = [];    % Corresponding k values
        k = 0;
        while true
            T_k = ((atan(d_delta_i(2)/d_delta_i(1)) + offset*pi) - w + k*pi)/(n + kappa*(eta*P + Q));
            if T_k > tf
                break;
            end
            T_opt_vals(end+1) = T_k;
            k_vals(end+1) = k;
            k = k + 1;
        end
        num_options = length(T_opt_vals);
        dv_vec = zeros(num_options, 3);
        for j=1:num_options
            if rem(k_vals(j) + offset, 2) == 0
                dv_vec(j, :) = [0, 0, i_dv_min]; % Normal is optimal
            else
                dv_vec(j, :) = [0, 0, -i_dv_min]; % Normal is optimal
            end
        end

        norms = vecnorm(dv_vec, 2, 2);    
        [~, idx] = min(norms);     
        delta_v_vals_inc = dv_vec(idx, :);
        delta_v_times_inc = T_opt(idx);
        
        fprintf('Appending inclination change burn to delta_v_vals...\n');
        delta_v_vals = [delta_v_vals; delta_v_vals_inc];
        delta_v_times = [delta_v_times; delta_v_times_inc];
    end

end
