function SimRunner_par2(RUNF, EXPNAME, RESUME_SIMULATION)

% SIMRUNNER_PAR Run experiments in parallel.
%    SIMRUNNER_PAR(RUNF, EXPNAME, RESUME_SIMULATION) requires that the 
%    (strictly named) main file 'runModel.m' is on the path and, like 
%    'test_main.m', accepts a structure of inputs 'inputs'.  
% 
% Mandatory Inputs
%    RUNF is the runfile containing the structured array CONSTANT and 
%    vectors of experimental conditions contained in structured array 
%    VARIABLE.  
%
% Optional Inputs
%    EXPNAME is the output .mat file name to save. Default: .mat version of
%    RUNF
%    RESUME_SIMULATION: whether to pick up where previous run crashed.
%    Determines pick-up point from .mat file. If no mat file exists, will
%    run from beginning.
%
%    Saves the same structure as in SIMRUNNER. Parfor run in batches of
%    replicates of each parameter case. So intermediate results can be
%    saved if parpool crash
%
%       Example: run the test runfile in parallel
%          SimRunner_par('test_runfile.txt')
%
%    For more information on setup see the helpfile for SIMRUNNER.
%
%    Notes:
%       - requires Distributed Computing Toolbox.
%
% See also SIMRUNNER MAKE_EXPS_TABLE SMARTTIME PREDICT_TIME

% Author: SA, 2 Oct 2017 (Monash)
% Modified: TM, 26 Jul 2024 (Curtin)

    %% Execution Setup
time_start = datetime();
print_progress = true; % turn to false for silent simulation
print_period_short = minutes(10); % print after every parameter case until this period is reached
print_freq_short = seconds(1);
print_freq_long = hours(1); % after reaching this value, reduce print frequency to once every longterm period
print_freq_medium = print_period_short;
print_period_current = "short"; % tracks whether currently in short, medium, or long term print period
print_freq_current = print_freq_short;
   
    %% Checks
if ~exist('EXPNAME', 'var') || isempty(EXPNAME)
    EXPNAME = strrep(RUNF, '.txt', '.mat');
end
if ~exist('RESUME_SIMULATION', 'var') || isempty(RESUME_SIMULATION') || ~exist(EXPNAME, 'file')
    RESUME_SIMULATION = false;
end

% .. check license
DO_PARALLEL = true;
if license('checkout', 'distrib_computing_toolbox') == 0
    fprintf(' !! Did not find distributed computing toolbox installed. Will use serial processing instead.\n');
    DO_PARALLEL = false;
end

    %% Simulation Setup
% -- Ingest constants and variables
[T, grouping_vars] = make_exps_table2(RUNF);

% -- Setup
% Notes:
%  . to satisify variable transparency, we pre-create the summary structure
%  . in addition, we must enforce the main proc, 'runModel()' since we cannot use `eval()` since it breaks transparency
crash_limit = 5;
n_ex = height(T);
n_rep = T{1, 'n_replicates'}; % assuming table sorted by RNG seed last
results_mode = string(T{1,'results_mode'});
network_plot_times = T{1, 'network_plot_times'};
rolling_average_window = T{1, 'rolling_average_window'};
% run_func = @(inputs) runModel(inputs);

n_cases = n_ex / n_rep;
a = struct('inputs', [], 'results', [], 'start', [], 'duration', []);
summary_case_empty = repmat(a, n_rep, 1);
if RESUME_SIMULATION % assuming M contains only completed simulations
    load(EXPNAME, 'M', 'S', 'inputs_collated');
    case_start = height(M) + 1; % each row in M is a aparameter case
else
    case_start = 1;
    if results_mode == "minimal" % saving only averaged results for plotting
        M = table();
        S = table();
    else % saving full timeseries
        summary = repmat(a, n_ex, 1);
    end
end
if print_progress
    fprintf('Setup finished: %s\n', datetime() - time_start);
end

%% -- Run (parallel)
if print_progress
    fprintf('Starting short-term progress printing ...\n');
    time_last_print = datetime();
end
for i_case = case_start :n_cases % running in batches of replicates of each parameter case
    ex_range = (i_case - 1) * n_rep + 1 : i_case * n_rep; % range of case replicate indices in full set
    T_case = T(ex_range,:);
    summary_case = summary_case_empty; % reset summary_case
        
    crash_count = 0; %0 % number of times parfor has crashed    
    while crash_count >= 0 && crash_count < crash_limit % retry upon error, but give up if re-occuring
        try % run parfor
            if DO_PARALLEL
                parfor ex = 1:n_rep
                    rng(1)              % for reproducibility
                    % .. get inputs for this run
                    inputs = table2struct(T_case(ex,:));
                    % .. save inputs, evaluate model
                    summary_case(ex).inputs  = inputs;
                    summary_case(ex).start   = datetime;                
                    summary_case(ex).results = runModel(inputs);
                    summary_case(ex).duration    = datetime() - summary_case(ex).start;
                end
            else
                for ex = 1:n_rep
                    rng(1)              % for reproducibility
                    % .. get inputs for this run
                    inputs = table2struct(T_case(ex,:));
                    % .. save inputs, evaluate model
                    summary_case(ex).inputs  = inputs;
                    summary_case(ex).start   = datetime;                
                    summary_case(ex).results = runModel(inputs);
                    summary_case(ex).duration    = datetime() - summary_case(ex).start;
                end
            end
        catch ME
            crash_count = crash_count + 1; 
            if print_progress
                fprintf('Crash #%d: %s\n', crash_count, ME.identifier);
            end
            if crash_count == crash_limit % give up on trying to get it to run
                fprintf('Unable to run without error.\n');
                rethrow(ME)
            end
            if strcmp(ME.identifier, 'parallel:lang:parfor:SessionShutDown') % if error is due to parpool crash
                delete(gcp('nocreate')); % shut down parpool
                pause(2^crash_count); % gap between retries widens, in case crash is due to current RAM/CPU load
                parpool(); % restart parpool
                continue; % retry parfor           
            else
                rethrow(ME) % if not expected parfor error, terminate MATLAB as per exception
            end
        end
        crash_count = -1; % indicates successful run
    end        
    if results_mode == "minimal"
        if crash_count < crash_limit
            [M_case, S_case, inputs_collated] = average_results_endstate(summary_case, grouping_vars);
            M = [M ; M_case]; % appending new rows
            S = [S ; S_case];
            % inputs_collated is corrected at end
        else
            M{end+1, 'n_replicates'} = n_rep; % adding dummy row for skipped case
            S{end+1, 1} = 0; % assuming first column is numeric entry
        end
    else
        summary(ex_range) = summary_case;
    end
    if print_progress
        if print_period_current == "short" && datetime() - time_start > print_period_short
            fprintf('Switching to medium-term printing ...\n');
            print_period_current = "medium";
            print_freq_current = print_freq_medium;            
        end
        if print_period_current == "medium" && datetime() - time_start > print_freq_long
            fprintf('Switching to long-term printing ...\n');
            print_period_current = "long";
            print_freq_current = print_freq_long;            
        end
        if datetime() - time_last_print > print_freq_current            
            fprintf('\tParameter Case #%d/%d (%.0f%%) completed: %s\n', i_case, n_cases, i_case/n_cases*100 ... 
                , datetime());
            time_last_print = datetime();
        end
    end
end
if print_progress
    fprintf('Runs finished: %s\n', datetime() - time_start);
end

    %% .. finish up, save
if results_mode == "minimal"
    inputs_collated.rng_seed = unique(T.rng_seed);
    save(EXPNAME,'M', 'S', 'inputs_collated');
elseif results_mode == "timeseries_plot"
    [inputs, res, T, G, NodeActions] = average_results_timeseries(summary, network_plot_times, rolling_average_window);
    save(EXPNAME, 'inputs', 'res', 'T', 'G', 'NodeActions');
else
    save(EXPNAME,'summary','-v7.3'); % for large timeseries
end

if print_progress
    fprintf('Program finished: %s\n', datetime - time_start);
end
