%% AA279D Simulations for PSET 4
% Tycho Bogdanowitsch
% Anshuk Chigullapalli 

%%%%% INITIAL RELATIVE ORBITAL ELEMENTS %%%%%%% 
% From PS3
d_a_SV2_init_1 = 0; % m
d_lambda_SV2_init_1 = 0; % m % -124
d_e_x_SV2_init_1 = 0; % m
d_e_y_SV2_init_1 = 100; % m
d_i_x_SV2_init_1 = 0; % m % 79
d_i_y_SV2_init_1 = 1000; % m % 1005

d_a_SV3_init_1 = 0; % m
d_lambda_SV3_init_1 = 0; % m % -793
d_e_x_SV3_init_1 = 0; % m
d_e_y_SV3_init_1 = 200; % m
d_i_x_SV3_init_1 = 0; % m
d_i_y_SV3_init_1 = 800; % m % 827

% Given in PS4
d_a_SV2_init_2 = 0; % m
d_lambda_SV2_init_2 = 100; % m % -124
d_e_x_SV2_init_2 = 50; % m
d_e_y_SV2_init_2 = 100; % m
d_i_x_SV2_init_2 = 30; % m % 79
d_i_y_SV2_init_2 = 200; % m % 1005

% Applying similar elements for SV3
d_a_SV3_init_2 = 0; % m
d_lambda_SV3_init_2 = -100; % m % -793
d_e_x_SV3_init_2 = 50; % m
d_e_y_SV3_init_2 = 100; % m
d_i_x_SV3_init_2 = 30; % m
d_i_y_SV3_init_2 = 200; % m % 827



% Define deputy initial conditions for ROE2ECI
deputy_inits = {
    'SV2_1', d_a_SV2_init_1, d_lambda_SV2_init_1, d_e_x_SV2_init_1, d_e_y_SV2_init_1, d_i_x_SV2_init_1, d_i_y_SV2_init_1;
    'SV3_1', d_a_SV3_init_1, d_lambda_SV3_init_1, d_e_x_SV3_init_1, d_e_y_SV3_init_1, d_i_x_SV3_init_1, d_i_y_SV3_init_1;
    'SV2_2', d_a_SV2_init_2, d_lambda_SV2_init_2, d_e_x_SV2_init_2, d_e_y_SV2_init_2, d_i_x_SV2_init_2, d_i_y_SV2_init_2;
    'SV3_2', d_a_SV3_init_2, d_lambda_SV3_init_2, d_e_x_SV3_init_2, d_e_y_SV3_init_2, d_i_x_SV3_init_2, d_i_y_SV3_init_2;
    'SV2_3', d_a_SV2_init_2, d_lambda_SV2_init_2, d_e_x_SV2_init_2, d_e_y_SV2_init_2, 0,                d_i_y_SV2_init_2;
    'SV3_3', d_a_SV3_init_2, d_lambda_SV3_init_2, d_e_x_SV3_init_2, d_e_y_SV3_init_2, 0,                d_i_y_SV3_init_2;
};

% Preallocate structs
state_abs = struct();
sim_results = struct();

% Chief initial state
state_abs.SV1 = [r_SV1_ECI_init; v_SV1_ECI_init];

% Loop through deputy initializations
for k = 1:size(deputy_inits,1)
    name = deputy_inits{k,1};
    d_a = deputy_inits{k,2};
    d_lambda = deputy_inits{k,3};
    d_e_x = deputy_inits{k,4};
    d_e_y = deputy_inits{k,5};
    d_i_x = deputy_inits{k,6};
    d_i_y = deputy_inits{k,7};

    % Convert to absolute ECI position/velocity
    [r_init, v_init] = ROE2ECI(a_SV1_init, ex_SV1_init, ey_SV1_init, i_SV1_init, RAAN_SV1_init, u_SV1_init, ...
                               d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y);

    state_abs.(name) = [r_init; v_init];
end

%%%%% RUN SIMULATIONS %%%%%

% List of spacecraft to simulate
sc_names = fieldnames(state_abs);

% Absolute dynamics simulations
for k = 1:length(sc_names)
    sc = sc_names{k}; % 'SV1', 'SV2_1', 'SV3_1', 'SV2_2', 'SV3_2'
    x0 = state_abs.(sc);

    % Simulate without J2
    [t_no_j2, state_no_j2] = rk4_eom_ECI(tstart:tint:tend, x0, false);
    sim_results.(sc).no_j2.t = t_no_j2;
    sim_results.(sc).no_j2.r = state_no_j2(:,1:3);
    sim_results.(sc).no_j2.v = state_no_j2(:,4:6);

    % Simulate with J2
    [t_with_j2, state_with_j2] = rk4_eom_ECI(tstart:tint:tend, x0, true);
    sim_results.(sc).with_j2.t = t_with_j2;
    sim_results.(sc).with_j2.r = state_with_j2(:,1:3);
    sim_results.(sc).with_j2.v = state_with_j2(:,4:6);
end

%%%%%%% CALCULATE OSC AND MEAN OE %%%%%%%%%%%

% Define cases
roe_cases = {
    'initial1', 'no_j2', 'SV2_1', 'SV3_1';
    'initial1', 'with_j2', 'SV2_1', 'SV3_1';
    'initial2', 'no_j2', 'SV2_2', 'SV3_2';
    'initial2', 'with_j2', 'SV2_2', 'SV3_2';
    'initial3', 'no_j2', 'SV2_3', 'SV3_3';
    'initial3', 'with_j2', 'SV2_3', 'SV3_3';
};

roe_results = struct();

% Loop over each case
for idx = 1:size(roe_cases,1)
    init_case = roe_cases{idx,1};  % 'initial1' or 'initial2'
    j2_case = roe_cases{idx,2};    % 'no_j2', 'with_j2'
    
    SV2_name = roe_cases{idx,3}; % 'SV2_1' or 'SV2_2'
    SV3_name = roe_cases{idx,4}; % 'SV3_1' or 'SV3_2'

    chief_name = 'SV1'; % Chief is always SV1

    % Grab chief position/velocity histories
    r_o = sim_results.(chief_name).(j2_case).r;
    v_o = sim_results.(chief_name).(j2_case).v;

    % List of deputies to run (both SV2 and SV3)
    deputies = {SV2_name, SV3_name};

    for deputy_idx = 1:2
        deputy_name = deputies{deputy_idx};

        % Grab deputy position/velocity histories
        r_t = sim_results.(deputy_name).(j2_case).r;
        v_t = sim_results.(deputy_name).(j2_case).v;

        % J2 flag for compute function
        if strcmp(j2_case, 'with_j2')
            J2_flag = 1;
        else
            J2_flag = 0;
        end

        [d_a_osc, d_lambda_osc, d_e_x_osc, d_e_y_osc, d_i_x_osc, d_i_y_osc, ...
         d_a_mean, d_lambda_mean, d_e_x_mean, d_e_y_mean, d_i_x_mean, d_i_y_mean, ...
         a_osc, e_x_osc, e_y_osc, i_osc, RAAN_osc, u_osc, ...
         a_mean, e_x_mean, e_y_mean, i_mean, RAAN_mean, u_mean] = ...
            compute_OE_ROE_mean_osc(r_o, v_o, r_t, v_t, J2_flag);

        
        [rho_RTN, rho_RTN_dot] = ECI2RTN_rel(r_o, v_o, r_t, v_t);
        
        % Save into the struct under a field for the deputy
        roe_results.(init_case).(j2_case).(deputy_name) = struct( ...
            'd_a_osc', d_a_osc, ...
            'd_lambda_osc', d_lambda_osc, ...
            'd_e_x_osc', d_e_x_osc, ...
            'd_e_y_osc', d_e_y_osc, ...
            'd_i_x_osc', d_i_x_osc, ...
            'd_i_y_osc', d_i_y_osc, ...
            'd_a_mean', d_a_mean, ...
            'd_lambda_mean', d_lambda_mean, ...
            'd_e_x_mean', d_e_x_mean, ...
            'd_e_y_mean', d_e_y_mean, ...
            'd_i_x_mean', d_i_x_mean, ...
            'd_i_y_mean', d_i_y_mean, ...
            'a_osc', a_osc, ...
            'e_x_osc', e_x_osc, ...
            'e_y_osc', e_y_osc, ...
            'i_osc', i_osc, ...
            'RAAN_osc', RAAN_osc, ...
            'u_osc', u_osc, ...
            'a_mean', a_mean, ...
            'e_x_mean', e_x_mean, ...
            'e_y_mean', e_y_mean, ...
            'i_mean', i_mean, ...
            'RAAN_mean', RAAN_mean, ...
            'u_mean', u_mean, ...
            'rho_RTN', rho_RTN, ...
            'rho_RTN_dot', rho_RTN_dot...
        );
    end
end

roe_analytical_results = struct();
deputy_mapping = struct( ...
    'SV2_1', 'initial1', ...
    'SV3_1', 'initial1', ...
    'SV2_2', 'initial2', ...
    'SV3_2', 'initial2', ...
    'SV2_3', 'initial3', ...
    'SV3_3', 'initial3' ...
);
% Loop through deputy initializations
for k = 1:size(deputy_inits,1)
    deputy_name = deputy_inits{k,1};
    d_a = deputy_inits{k,2};
    d_lambda = deputy_inits{k,3};
    d_e_x = deputy_inits{k,4};
    d_e_y = deputy_inits{k,5};
    d_i_x = deputy_inits{k,6};
    d_i_y = deputy_inits{k,7};

    condition = deputy_mapping.(deputy_name);

    deputy_ic_given = [d_a, d_lambda, d_e_x, d_e_y, d_i_x, d_i_y];
    % TODO: FIGURE OUT WHAT IS HAPPENING HERE
    deputy_output = roe_results.(condition).with_j2.(deputy_name);
    deputy_ic_mean = [deputy_output.d_a_mean(1), deputy_output.d_lambda_mean(1), deputy_output.d_e_x_mean(1), deputy_output.d_e_y_mean(1), deputy_output.d_i_x_mean(1), deputy_output.d_i_y_mean(1)];

    chief_name = 'SV1'; % Chief is always SV1
    t_series = sim_results.(chief_name).('no_j2').t; % Pulling a random t vector
    SV1_oe_init = [a_SV1_init, e_SV1_init, i_SV1_init, RAAN_SV1_init, w_SV1_init];
    roe_output_given_ic = roe_stm_j2(t_series, deputy_ic_given, SV1_oe_init);
    roe_output_mean_ic = roe_stm_j2(t_series, deputy_ic_mean, SV1_oe_init);

    roe_analytical_results.(condition).(deputy_name) = struct('roe_analytical_j2_given_ic', roe_output_given_ic, ...
                                                              'roe_analytical_j2_mean_ic', roe_output_mean_ic); % Using the initial conditions as given (0 delta a) and mean ICs (non-zero delta a)
end

plot_OE_ROE_mean_osc(roe_results, tstart, tint, tend, t_orbit)
plot_ROE_compare_analytical(roe_results, roe_analytical_results, tstart, tint, tend, t_orbit);

%%%%%%%%%%%%% CONVERT TO RTN %%%%%%%%%%%%%%%%
% Convert the absolute ECI positions of the chief and deputies and convert to RTN.
%[rho_SV2_RTN, rho_SV2_RTN_dot] = ECI2RTN_rel(r_SV1_ECI_no_j2, v_SV1_ECI_no_j2, r_SV2_ECI_no_j2, v_SV2_ECI_no_j2);
%[rho_SV3_RTN, rho_SV3_RTN_dot] = ECI2RTN_rel(r_SV1_ECI_no_j2, v_SV1_ECI_no_j2, r_SV3_ECI_no_j2, v_SV3_ECI_no_j2);