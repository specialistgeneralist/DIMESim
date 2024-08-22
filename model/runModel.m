function OUT = runModel(IN)

%% Run main DIME simulation

% -- control rng
rng(IN.rng_seed);

% -- Set up global state vars
global U D I M E A C xh x
global G_Outcome_D G_Outcome_I G_Outcome_M G_Outcome_E
global G_Type_D G_Type_I G_Type_M G_Type_E
global G_Interaction_D G_Interaction_I G_Interaction_M G_Interaction_E

% -- Set up additional global state vars as required for variant
if ~contains(IN.model_variant, 'GeneralBroadcast')
    global B_after_IR 
end
if contains(IN.model_variant, 'CollectiveReinterpretation')
    global B_after_IR_succ_bool B_after_CR S
    global adjMatrix pMatrix          % pMatrix(i,j) = 1/deg(i) if i-j is an edge
end

% -- Initialise model (from static file, and params)
IN.studyParams = get_studyParams(IN.study_params_file);   % -- read in study params
InitialiseModel(IN);                                      % -- init all global vars

% .. timing vars
last_active_time = zeros(IN.n, 1); % xi in manuscript
t_stop = [];

%% Enter main loop
stopping_cond_met = false;
for tt = 1 : IN.T

    % // Generate broadcast
    U(tt + 1) = broadcast(IN.p);
    Bf = U(tt + 1);

    % -- apply interpretation based on informational variant
    if ~contains(IN.model_variant, 'GeneralBroadcast')                   % -- do for everything other than General Broadcast
        % Perceived broadcast after individual reinterpretation
        B_after_IR(:, tt + 1) = perceived_success_after_IR_vec(U(tt + 1), D(:, tt), IN.F);
        Bf = B_after_IR(:, tt + 1);
    end
    if contains(IN.model_variant, 'CollectiveReinterpretation')
        B_after_CR_old = zeros(size(B_after_CR));
        %tt


        % Implement R rounds of the CR step
        B_after_CR_old(:, tt + 1) = B_after_IR(:, tt + 1);
        for jj = 1:IN.R

            % Reinitialise variables
            B_after_IR_succ_bool(:,tt + 1) = zeros(IN.n,1);
            S(:,tt + 1) = zeros(IN.n, 1);
            B_after_CR(:, tt + 1) = zeros(IN.n, 1);


            % Calculate  number of successes in neighbourbood
            B_after_IR_succ_bool(:, tt + 1) = (1/2) .* (ones(IN.n,1) - B_after_CR_old(:, tt + 1));
            S(:, tt + 1) = pMatrix * B_after_IR_succ_bool(:, tt + 1);

            % Perceived broadcast after collective reinterpretation
            B_after_CR(:, tt + 1) = perceived_success_after_CR_vec(S(:, tt + 1), B_after_IR_succ_bool(:, tt + 1), IN.nu);
            
            %
            B_after_CR_old(:, tt + 1) = B_after_CR(:, tt + 1);
            
            %jj
            %B_after_CR_old(:, tt + 1)
        end
        Bf = B_after_CR(:, tt + 1);

    end

    % // Evolve DIME variables, action intention and current action focus values, using array algebra (i.e. simultaneously)

    % -- Choose x for regression
    switch lower(IN.x_regression_type)
        case 'action'
            x_reg = x(:, tt);
        case 'orientation'
            x_reg = C(:, tt) .* xh(:, tt);
        case 'last active action'
            x_reg = xh(:, tt);        
    end    

    % -- Evolve DIME variables    
    D(:, tt + 1) = p_sat(D(:, tt) + Bf .* G_Outcome_D + x_reg .* G_Type_D + x_reg .* Bf .* G_Interaction_D + IN.DIME_update_temp * unifrnd(-1,1, [IN.n, 1]));
    I(:, tt + 1) = p_sat(I(:, tt) + Bf .* G_Outcome_I + x_reg .* G_Type_I + x_reg .* Bf .* G_Interaction_I + IN.DIME_update_temp * unifrnd(-1,1, [IN.n, 1]));
    M(:, tt + 1) = p_sat(M(:, tt) + Bf .* G_Outcome_M + x_reg .* G_Type_M + x_reg .* Bf .* G_Interaction_M + IN.DIME_update_temp * unifrnd(-1,1, [IN.n, 1]));
    E(:, tt + 1) = p_sat(E(:, tt) + Bf .* G_Outcome_E + x_reg .* G_Type_E + x_reg .* Bf .* G_Interaction_E + IN.DIME_update_temp * unifrnd(-1,1, [IN.n, 1]));
    
    % -- Evolve action intention and current action focus values
    A(:, tt + 1) = action_intention(D(:, tt + 1), I(:, tt + 1), M(:, tt + 1), E(:, tt + 1));
    C(:, tt + 1) = action_focus(I(:, tt + 1), M(:, tt + 1), E(:, tt + 1));

        
    % -- Evolve historical action
    xh(:, tt + 1) = xh(:, tt); % initially assume last action is same from previous time step    
    active_previous_time = find(x(:, tt) ~= 0); % agents who were active in the previous time step
    xh(active_previous_time, tt + 1) = x(active_previous_time, tt); % update last action for these 

    % -- Evolve action
    x(:, tt + 1) = action(A(:, tt + 1), C(:, tt + 1), xh(:, tt + 1));

    % Check for stopping conditions, based on D I M E vars entering epsilon neighbourhood of last w steps., if met, end run
    % -- note we only check every stopping_check_interval steps to save compute
    if IN.do_stopping && tt >= IN.stopping_window_steps + IN.rolling_average_window - 1 ...
        && mod(tt, IN.stopping_check_interval) == 0
        stopping_cond_met = check_stopping_conditions(tt, IN);
    end
    if stopping_cond_met
        t_stop = tt;
        break
    end

end

% -- Gather results
% .. if we are in 'full_run' mode, take full time-series up till stopping
% .. if in 'end_state_only' mode, take last `stopping_window_steps` steps only
switch IN.results_mode
    case {'full_run', 'timeseries_plot'}
        t_record_cols = 1:tt+1;
    case {'end_state_only', 'minimal'}
        t_record_cols = tt+1 - IN.stopping_window_steps + 1:tt+1;    
end
if ~stopping_cond_met
    t_stop = tt;
end

% .. broadcasts
OUT.U = single(U(t_record_cols));
if ~contains(IN.model_variant, 'GeneralBroadcast')
    OUT.B_after_IR = single(B_after_IR(:, t_record_cols));
end
if contains(IN.model_variant, 'CollectiveReinterpretation')
    OUT.B_after_CR = single(B_after_CR(:, t_record_cols));
end
% .. DIME vars
OUT.D = single(D(:, t_record_cols));
OUT.I = single(I(:, t_record_cols));
OUT.M = single(M(:, t_record_cols));
OUT.E = single(E(:, t_record_cols));
% .. action vars
OUT.A = single(A(:, t_record_cols));
OUT.C = single(C(:, t_record_cols));
OUT.xh = single(xh(:, t_record_cols));
OUT.x = single(x(:, t_record_cols));
% .. other
OUT.t_stop = t_stop;
OUT.stopping_cond_met = stopping_cond_met;
OUT.adjMatrix = adjMatrix;

end
