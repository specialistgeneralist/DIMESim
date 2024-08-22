function simrunner(MAINFXN, RUNF)

% SIMRUNNER Run full-factorial numerical experiments, in series.
%    SIMRUNNER(MAINFXN, RUNF) Runs MAINFXN over the constants and 
%    full-factorial combinations of variables found in RUNF, creating
%    an experiment EXPNAME from the stem of RUNF, and saving the 
%    results as a structured array SUMMARY to the output file EXPNAME.mat.
%
%    During run-time, progress is provided to the workspace, and, with 
%    human-readable estimates of the overall finish time. Finish time is 
%    computed after a few runs using ensemble logistic regression trees 
%    (FITRENSEMBLE) to make allowances for the likely run-time of certain 
%    future combinations of variables. If the statistics_toolbox is not 
%    installed, the prediction will fall back on a simple average of runs 
%    already completed.
%
%       Example: run the test program and runfile
%          simrunner('test_main', 'test_runfile.txt')
%
%    For a parallel version of SIMRUNNER, see SIMRUNNER_PAR.
%
%    Setup:
%       main program
%          The main program file (MAINFXN) should take a single structured 
%          array as input. This input array will be created at run-time by 
%          SIMRUNNER, as a combination of constants and a unique parameter set 
%          from the variables given. Params can be easily unpacked from the 
%          input array within the main file. Note: SIMRUNNER can handle both 
%          vectors and cell-array of strings as variable definitions.
%
%       runfile
%          The runfile (RUNF) should be a script, containing the definitions 
%          of exactly two structured arrays, CONSTANT and VARIABLE. For 
%          example:
%              constant.report_flag = 0;
%              constant.temp        = 25.5;
%              variable(1).name     = 'a';
%              variable(1).values   = {'red' 'green' 'blue'};
%              variable(2).name     = 'b';
%              variable(2).values   = [2 4 8 16];
%          This file describes two constants, and two variables named 'a' and 
%          'b', with 3 and 4 distinct values to be tested respectively. As a 
%          full-factorial experiment, the 12 unique combinations will be run 
%          (3x4), resulting in an output structure with 12 entries.
%
%          Note: if the constant structure is more than 1 depth, estimates of 
%          goodness of fit for time remaining in experiment will not be 
%          possible since full predictor values will not be transparent to the 
%          ensemble tree fitting algorithm.
%
% See also SIMRUNNER MAKE_EXPS_TABLE SMARTTIME PREDICT_TIME

% Author: SA, 2 Oct 2017 (Monash)

more off

% .. Ingest constants and variables
T = make_exps_table(RUNF);

% .. build experiment name
EXPNAME = strrep(RUNF,'.txt','');

% .. Setup
n_ex = height(T);
runtime = [];

% .. Run
for ex = 1:n_ex
    rng(1);         % for reproducibility
    % .. get inputs for this run
    inputs = table2struct(T(ex,:));
    % .. update info
    clc, type message.txt
    fprintf('  >> Running ''%s'' > experiment ''%s''\n', MAINFXN, EXPNAME);
    fprintf('  >> Simulation continuing: %d of %d\n',ex,n_ex);
	% report on progress (time, and visually)
    if n_ex > 1
        if ex > 2
            [eta_num,fit] = predict_time(ex-1,runtime,T);
            eta_rec(ex) = eta_num;
			eta = datestr(eta_num);
			smart_eta = smarttime(eta);
            fprintf('  >> Estimated final finish time (fit):\n          ...  %s (%s)\n ',smart_eta,fit);
        else
            fprintf('  >> Estimated final finish time: <more data needed>\n');
		end
        P = visprog(ex,n_ex,50,1);
		disp(P);
    else
        fprintf('  >> Single experiment running (no prediction available)');
	end
    % .. save inputs, evaluate model
    summary(ex).inputs = inputs;
    summary(ex).start = datestr(now);
    eval(['summary(ex).results=' MAINFXN '(inputs);']);
    summary(ex).stop = datestr(now);
    runtime(ex,:) = datenum(summary(ex).stop) - datenum(summary(ex).start);
	% finish up
    save([EXPNAME],'summary');
end
clc, type endmessage.txt
clear all
