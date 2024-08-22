function simrunner_par(RUNF)

% SIMRUNNER_PAR Run experiments in parallel.
%    SIMRUNNER_PAR(RUNF, EXPNAME) requires that the (strictly named) main file 
%    'runModel.m' is on the path and, like 'test_main.m', accepts a structure 
%    of inputs 'inputs'.  RUNF is the runfile containing the structured array 
%    CONSTANT and vectors of experimental conditions contained in structured 
%    array VARIABLE.  EXPNAME is the output .mat file name to save.
%
%    Saves the same structure as in SIMRUNNER, however, due to the parallel 
%    nature of this version, only saves the outcome once all parallel jobs are 
%    finished.
%
%       Example: run the test runfile in parallel
%          SimRunner_par('test_runfile.txt')
%
%    For more information on setup see the helpfile for SIMRUNNER.
%
%    Notes:
%       - no information on progress is provided in this version.
%       - requires Distributed Computing Toolbox.
%
% See also SIMRUNNER MAKE_EXPS_TABLE SMARTTIME PREDICT_TIME

% Author: SA, 2 Oct 2017 (Monash)

% TODO
%  SA: Progress meter. This is not trivial in parallel.

% .. check license
if license('checkout', 'distrib_computing_toolbox') == 0
    error('Did not find distributed computing toolbox installed.')
end

% -- Ingest constants and variables
T = make_exps_table(RUNF);

% -- Setup
% Notes:
%  . to satisify variable transparency, we pre-create the summary structure
%  . in addition, we must enforce the main proc, 'runModel()' since we cannot use `eval()` since it breaks transparency
n_ex = height(T);
runtime = [];
a = struct('inputs', [], 'results', [], 'start', [], 'stop', []);
summary = repmat(a, n_ex, 1);
EXPNAME = strrep(RUNF, '.txt', '.mat');

% -- Run (parallel)
parfor ex = 1:n_ex
    rng(1)              % for reproducibility
    % .. get inputs for this run
    inputs = table2struct(T(ex,:));
    % .. save inputs, evaluate model
    summary(ex).inputs  = inputs;
    summary(ex).start   = datestr(now);
    summary(ex).results = runModel(inputs);
    summary(ex).stop    = datestr(now);
end
% .. finish up, save
save([EXPNAME],'summary');
