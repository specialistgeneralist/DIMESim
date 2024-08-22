function T = make_exps_table(RUNF)

% MAKE_EXPS_TABLE Generate full factorial variable table.
%    T = MAKE_EXPS_TABLE(RUNF) creates a table T for full factorial 
%    exploration of a model. RUNF is the name of a file (e.g. 'runfile.txt') 
%    on the path that provides two structs, CONSTANT, the constants of the 
%    experiments, and VARIABLE, the variables to explore.
%
%    Note: VARIABLE is optional. If no VARIABLE structure is provided, a 
%    single row table will be returned, the single experiment to be run with 
%    each constant as cols.
%
%    VARIABLE should have the following structure:
%       variable(1).name = 'a';
%       variable(1).values = {'red' 'green' 'blue'};
%       variable(2).name = 'b';
%       variable(2).values = [2 4 8 16];
%       ...
%    T will provide a table of rows PROD(num_opts per variable) and one col 
%    per variable, named by variable names, and one col per constant.
%
%    Examples:
%        T = make_exps_table('test_runfile.txt')
%
%    See also NDGRID

% Author: SA, 2 Oct 2017 (Monash)

% TODO
%  SA: add support for semi-factorial design

% .. check we can find the RUNF
if exist(RUNF,'file') == 0
    error('Specified runfile not found.');
end

% .. read runfile
fid = fopen(RUNF,'r');
Rnf = char([fread(fid)]');     % --> 'constant' 'variable'
fclose('all');
eval([Rnf']);

% .. check inputs exists
if exist('constant', 'var') == 0
    error('Found runfile, but doesn''t seem to contain structured array variable ''constant''.');
end

% .. assume only constants, unless variables found
n = 1;
if exist('variable', 'var')
    % .. create full factorial experiment index table
    lvars = length(variable);
	for i = 1:lvars
		nvars(i) = size(variable(i).values,2);
	    if i == 1
	        rhargs = sprintf('1:%.0f', nvars(i));
	        lhargs = sprintf('[x1');
	    else
	        rhargs = sprintf('%s, 1:%.0f', rhargs, nvars(i));
            lhargs = sprintf('%s, x%.0f', lhargs, i);
	    end
	end
	ndgrid_s = sprintf('%s] = ndgrid(%s);', lhargs, rhargs);
    eval(ndgrid_s);     % create x1..xk

    % .. now turn into experimental table
    ix = 1:numel(x1);
    for i = 1:lvars
        tname = variable(i).name;
        eval(sprintf('ix_t = x%.0f(ix);', i));
        t1 = variable(i).values;
        t2 = [t1(ix_t)]';
        if i == 1
            T = table(t2);
        else
            T = [T table(t2)];
        end
        T.Properties.VariableNames{end} = tname;
    end
    T = unique(T, 'rows');      % .. hierarchy by var 1 .. k
    n = height(T);
end

% .. now add in constants
cnames = fieldnames(constant);
lconsts = numel(cnames);
for i = 1:lconsts
    this_c = cnames{i};
    if ischar(constant.(this_c))
        r = repmat({constant.(this_c)}, n, 1);
    else
        r = repmat(constant.(this_c), n, 1);
    end
    T.(this_c) = r;
end
% .. convert to table in case of no variables
if isstruct(T)
    T = struct2table(T);
end
