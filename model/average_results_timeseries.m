function [inputs, res, T, G, NodeActions] = average_results_timeseries(summary, NETWORK_PLOT_TIMES, ROLLING_AVERAGE_WINDOW)
%%% Plotting time series of DIME Attributes and Action Fractions
    tic;
    %% Constants and Inputs
% -- constant settings
% pc_var = @(x) (max(x) - min(x))/(mean(x) + (mean(x)==0)); % denominator modification ensures pc_var = 0 when mean is 0

% -- input parameters
% ALPHA = 0.05;           % -- confidence interval alpha
% NBOOT = 1000;           % -- number of bootstrap samples to use for confidence intervals
% NETWORK_PLOT_TIMES = [];    % -- times at which network should be plotted
% ROLLING_AVERAGE_WINDOW = []; % -- window for calculating rolling average


    %% Calculation Setup
n_npt = numel(NETWORK_PLOT_TIMES);
    
n = numel(summary);         % number of simulations, over both parameter sweeps and replicates
inputs = summary(1).inputs;
results = [summary(:).results];
stopping_window_steps = inputs.stopping_window_steps;
% stopping_check_interval = inputs.stopping_check_interval;
% stopping_epsilon = inputs.stopping_epsilon;
full_run_mode = ismember(inputs.results_mode, {'full_run', 'timeseries_plot'});
n_agents = inputs.n;
if isempty(ROLLING_AVERAGE_WINDOW)
    ROLLING_AVERAGE_WINDOW = stopping_window_steps;
end

% tFirstCheck = (floor(stopping_window_steps / stopping_check_interval) +1) * stopping_check_interval;
if full_run_mode
    T = inputs.T;
    t_converged_set = zeros(n,1); % will build up in next section
else    
    t_converged_set = sort(results([results(:).stopping_cond_met]).t_stop)';  % t_stop for only converged cases    
    T = max(t_converged_set); % find maximum time length needed to be plotted    
end
if isempty(NETWORK_PLOT_TIMES)
    NETWORK_PLOT_TIMES = [1 round(T/2) T];
end

    %% Extracting Results
res = repmat(struct('count', 0, 'ConHab_f', 0, 'RadHab_f', 0, 'Inv_f', 0 ...
    ,'InaConHab_f', 0, 'InaRadHab_f', 0, 'InaInv_f', 0 ...
    ,'D_f', 0, 'I_f', 0, 'M_f', 0, 'E_f', 0),T,1); % initializing struct array. Adding null field as array creation requires at least one existing field
for i = 1:n
    % -- get fraction of agents who are active (conventional, radical) and inactive
    Cs = results(i).C;        % -- subset of C: 1 if conventional, -1 if radical
    As = results(i).A;        % -- subset of A: 1 if active, 0 if inactive
    xhs = results(i).xh;
    xs =  results(i).x;
    D = results(i).D;
    I = results(i).I;
    M = results(i).M;
    E = results(i).E;       
    [x_con, x_rad, x_ina, xi_conhab, xi_coninv, xi_radhab, xi_radinv, xi_conlat, xi_radlat, x_conhab, x_coninv, x_radhab, x_radinv] = std_outputs_actions(As, Cs, xhs, xs);
    
    t_set = 1:results(i).t_stop;    
    for t = t_set
        res(t).count = res(t).count + 1;
        % res(t).Con_f = res(t).Con_f + x_con(t);
        % res(t).Rad_f = res(t).Rad_f + x_rad(t);
        % res(t).Ina_f = res(t).Ina_f + x_ina(t);
        % res(t).InaConLat_f = res(t).InaConLat_f + xi_conlat(t);
        % res(t).InaRadLat_f = res(t).InaRadLat_f + xi_radlat(t);
        % res(t).C_f = res(t).C_f + mean(Cs(:,t) == -1); % fraction of agents innovating; i.e. switching action focus this time step
        res(t).ConHab_f = res(t).ConHab_f + x_conhab(t);
        res(t).RadHab_f = res(t).RadHab_f + x_radhab(t);
        res(t).Inv_f = res(t).Inv_f + x_coninv(t) + x_radinv(t);
        res(t).InaConHab_f = res(t).InaConHab_f + xi_conhab(t);
        res(t).InaRadHab_f = res(t).InaRadHab_f + xi_radhab(t);
        res(t).InaInv_f = res(t).InaInv_f + xi_coninv(t) + xi_radinv(t);        
        res(t).D_f = res(t).D_f + mean(D(:,t)); % average over agents
        res(t).I_f = res(t).I_f + mean(I(:,t));
        res(t).M_f = res(t).M_f + mean(M(:,t));
        res(t).E_f = res(t).E_f + mean(E(:,t));        
                
        % if full_run_mode && t_converged_set(i) == 0 && t >= tFirstCheck
        %     t0 = max(1, t - stopping_window_steps + 1);
        %     if max([pc_var(mean(D(:, t0:t),1)) pc_var(mean(I(:, t0:t),1)) ...
        %             pc_var(mean(M(:, t0:t),1)) pc_var(mean(E(:, t0:t),1))]) ...
        %             < stopping_epsilon
        %         t_converged_set(i) = t;
        %     end
        % end
    end    
end
% if full_run_mode
%     t_converged_set = sort(t_converged_set(t_converged_set ~= 0));    
% end
% t_converged_unique = reshape(unique(t_converged_set), 1, []);
% n_converged_unique = numel(t_converged_unique);
% converged_cum_timeseries = (1:numel(t_converged_set))' / n * 100; % cumulative fraction of replicates converged at each t_stop
fprintf('Replicate calculations finished: %.2f\n', toc);

% T = size(res,1);
t_set = 1:T;
for t = t_set % averaging over replicates. Replicates can stop at different times, so individual counts are needed
    t0 = max(1, t - ROLLING_AVERAGE_WINDOW + 1);        
    % for stem_ID = ["Con" "Rad" "Ina" "InaConLat" "InaRadLat" "C" "D" "I" "M" "E"]        
    for stem_ID = ["ConHab" "RadHab" "Inv" "InaConHab" "InaRadHab" "InaInv" "D" "I" "M" "E"]        
        res(t).(stem_ID + "_f") = res(t).(stem_ID + "_f") / res(t).count;
        res(t).(stem_ID + "_ra") = mean([res(t0:t).(stem_ID + "_f")]);    % rolling average
        % res(t).(stem_ID + "_var") = pc_var([res(t0:t).(stem_ID + "_f")]); % variance
        % res(t).(stem_ID + "_ra_var") = pc_var([res(t0:t).(stem_ID + "_ra")]); % variance on rolling average
    end
    % res(t).max_ra_var = max([res(t).D_ra_var, res(t).I_ra_var res(t).M_ra_var, res(t).E_ra_var]);   
end
fprintf('Replicate averaging finished: %.2f\n', toc);

action_struct = struct('Type', "", 'A', [], 'C', [], 'xh', []);
i = 1;
action_struct(i) = struct('Type', "ConHab", 'A', 1, 'C', 1, 'xh', 1); i = i +1;
action_struct(i) = struct('Type', "Inv", 'A', 1, 'C', -1, 'xh', 1); i = i +1;
action_struct(i) = struct('Type', "Inv", 'A', 1, 'C', -1, 'xh', -1); i = i +1;
action_struct(i) = struct('Type', "RadHab", 'A', 1, 'C', 1, 'xh', -1); i = i +1;
action_struct(i) = struct('Type', "InaConHab", 'A', 0, 'C', 1, 'xh', 1); i = i +1;
action_struct(i) = struct('Type', "InaInv", 'A', 0, 'C', -1, 'xh', 1); i = i +1;
action_struct(i) = struct('Type', "InaInv", 'A', 0, 'C', -1, 'xh', -1); i = i +1;
action_struct(i) = struct('Type', "InaRadHab", 'A', 0, 'C', 1, 'xh', -1); i = i +1;

G = graph(results(1).adjMatrix); % arbitrarily choosing first replicate for network data
NodeActions = strings(n_npt, n_agents);
for i = 1: n_npt
    t = NETWORK_PLOT_TIMES(i);
    for j = 1:n_agents
        % NodeAction = '';
        k = find([action_struct(:).A] == results(1).A(j,t) ...
            & [action_struct(:).C] == results(1).C(j,t) ...
            & [action_struct(:).xh] == results(1).xh(j,t), 1);
        if isempty(k)
            NodeActions(i,j) = "InaConHab";
            disp('Warning: Invalid action value');
        else
            NodeActions(i,j) = action_struct(k).Type;
        end
        %{
        switch results(1).x(j,t)
            case 1
                NodeActions(i, j) = "Con";
            case -1
                NodeActions(i,j) = "Rad";
            case 0
                switch results(1).xh(j,t) * results(1).C(j,t)
                    case 1
                        NodeActions(i,j) = "InaConLat";
                    case -1
                        NodeActions(i,j) = "InaRadLat";
                    otherwise
                        NodeActions(i,j) = "InaConLat";
                        disp('Warning: Latent action of 0');
                end
            otherwise
               
        end        
        %}
    end
    %1 - (0 + 1 * degreeList / max(degreeList)) * (1 - defaultNodeColor);        
end