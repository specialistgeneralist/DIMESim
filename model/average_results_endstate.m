function [M, S, inputs_collated, inputs_scalar] = average_results_endstate(summary, GROUPING_VARS)

M = table();
S = table();
inputs_collated = struct('grouping_vars', GROUPING_VARS);
inputs_scalar = struct();
pc_var = @(x) (max(x) - min(x))/(mean(x) + (mean(x)==0)); % denominator modification ensures pc_var = 0 when mean is 0

inputs = [summary(:).inputs];
for var = string(fields(inputs)')
    switch class(inputs(1).(var))
        case {'char', 'string'}
            var_array = {inputs(:).(var)};
        case 'logical'
            var_array = boolean([inputs(:).(var)]);
        otherwise
            var_array = [inputs(:).(var)];
    end
    var_unique = unique(var_array);

    if numel(var_unique) == 1
        if isequal(class(var_unique), 'cell')
            inputs_collated.(var) = var_unique{1};
            inputs_scalar.(var) = var_unique{1};
        else
            inputs_collated.(var) = var_unique;
            inputs_scalar.(var) = var_unique;
        end
    else
        inputs_collated.(var) = var_unique;
        if ischar(var_unique)
            inputs_scalar.(var) = var_unique;
        end        
    end
end

n = numel(summary);         % number of simulations, over both parameter sweeps and replicates
res = repmat(struct([]),n,1); % initializing struct array. Adding null field as array creation requires at least one existing field

% -- get average of last N steps over all agents, for each variable, for each replicate
for i = 1:n
    stopping_window_steps = summary(i).inputs.stopping_window_steps;      % -- window to build mean state variables over
    stopping_epsilon = summary(i).inputs.stopping_epsilon; % -- if range(y)/mean(y) < epsilon, stop (for all DIME var simultaneously)
    for var = reshape(GROUPING_VARS,1,[])
        res(i).(var) = summary(i).inputs.(var);
    end

    % -- mean DIME state variables over LAST_N_STEPS
    res(i).D_avg = mean(summary(i).results.D(:, end-stopping_window_steps+1:end), 'all');
    res(i).I_avg = mean(summary(i).results.I(:, end-stopping_window_steps+1:end), 'all');
    res(i).M_avg = mean(summary(i).results.M(:, end-stopping_window_steps+1:end), 'all');
    res(i).E_avg = mean(summary(i).results.E(:, end-stopping_window_steps+1:end), 'all');

    % -- PC variance of mean DIME state variables over stopping_window_steps
    res(i).D_pcvar = pc_var(mean(summary(i).results.D(:, end-stopping_window_steps+1:end), 1));
    res(i).I_pcvar = pc_var(mean(summary(i).results.I(:, end-stopping_window_steps+1:end), 1));
    res(i).M_pcvar = pc_var(mean(summary(i).results.M(:, end-stopping_window_steps+1:end), 1));
    res(i).E_pcvar = pc_var(mean(summary(i).results.E(:, end-stopping_window_steps+1:end), 1)); 
    res(i).max_pcvar = max([res(i).D_pcvar, res(i).I_pcvar, res(i).M_pcvar, res(i).E_pcvar]);
    res(i).unconverged = res(i).max_pcvar > stopping_epsilon;

    % -- get fraction of agents who are active (conventional, radical) and inactive
    Cs = summary(i).results.C(:, end-stopping_window_steps+1:end);        % -- subset of C: 1 if conventional, -1 if radical
    As = summary(i).results.A(:, end-stopping_window_steps+1:end);        % -- subset of A: 1 if active, 0 if inactive
    xhs = summary(i).results.xh(:, end-stopping_window_steps+1:end);
    xs =  summary(i).results.x(:, end-stopping_window_steps+1:end);    
    [x_con, x_rad, x_ina, xi_conhab, xi_coninv, xi_radhab, xi_radinv, xi_conlat, xi_radlat, x_conhab, x_coninv, x_radhab, x_radinv] = std_outputs_actions(As, Cs, xhs, xs);

    res(i).Con_f = mean(x_con);
    res(i).Ina_f = mean(x_ina);
    res(i).Rad_f = mean(x_rad);
    res(i).ConHab_f = mean(x_conhab);
    res(i).ConInv_f = mean(x_coninv);
    res(i).RadHab_f = mean(x_radhab);
    res(i).RadInv_f = mean(x_radinv);
    res(i).Inv_f = res(i).ConInv_f + res(i).RadInv_f;        

    res(i).InaConHab_f = mean(xi_conhab);
    res(i).InaConInv_f = mean(xi_coninv);
    res(i).InaRadHab_f = mean(xi_radhab);
    res(i).InaRadInv_f = mean(xi_radinv);
    res(i).InaConLat_f = mean(xi_conlat);
    res(i).InaRadLat_f = mean(xi_radlat);
    res(i).InaInv_f = res(i).InaConInv_f + res(i).InaRadInv_f;
end

% -- averaging replicate experiments and computing conf intervals using bootstrapping
t = struct2table(res);
vars = t.Properties.VariableNames;
avg_vars = vars(contains(vars, {'_avg' '_f' '_pcvar' 'unconverged'}));
for i = 1:numel(avg_vars)
    this_var = avg_vars{i};

    % .. get mean :: produces pairs stacked vertically.
    mm = rowfun(@mean, t,...
        'InputVariables', this_var,...
        'GroupingVariables', GROUPING_VARS,...
        'OutputVariableNames', this_var);
    mm.GroupCount = [];

    % .. collect
    if i == 1
        M = mm;
    else
        M = join(M, mm, 'Keys', GROUPING_VARS, 'RightVariables', this_var);
    end
    
    if size(t,1) > 1   % compute confidence interval only when there are multiple values of the grouping variables
        this_var_ci = [this_var '_ci'];
        this_var_ci_low = [this_var '_ci_low'];
        this_var_ci_high = [this_var '_ci_high'];
        % .. get CI :: produces pairs stacked vertically.
        ss = rowfun(@(x) bootci(inputs_collated.nboot, {@mean, x}, 'Alpha', inputs_collated.alpha), t,...
            'InputVariables', this_var,...
            'GroupingVariables', GROUPING_VARS,...
            'OutputVariableNames', {this_var_ci});
        % .. unstack (bit fiddly)
        ss.GroupCount = [];
        ss_low  = ss(1:2:end,:);
        ss_high = ss(2:2:end,:);
        ss_low = renamevars(ss_low, this_var_ci, this_var_ci_low);
        ss_high = renamevars(ss_high, this_var_ci, this_var_ci_high);
        ss = join(ss_low, ss_high, 'Keys', GROUPING_VARS);
    
        % .. collect
        if i == 1
            S = ss;
        else
            S = join(S, ss, 'Keys', GROUPING_VARS, 'RightVariables', {this_var_ci_low this_var_ci_high});
        end
    end    
end

for var = string(fields(inputs_scalar)') % so that M contains all data relevant to results
    var_value = inputs_scalar.(var);
    if ischar(var_value)
        M.(var)(:) = {var_value};
    else
        M.(var)(:) = var_value;
    end
end
