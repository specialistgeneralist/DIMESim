function IN = read_constant_runfile(FNAME)

%% Read a single run-file for simple running of the model (constant only, no variables).
%
% NOTE: we assume no replicates or variables are being tested, just a single model run (one replicate).

fid = fopen(FNAME,'r');
Rnf = char([fread(fid)]');     % --> 'constant' 'variable'
eval([Rnf]');
IN = constant;