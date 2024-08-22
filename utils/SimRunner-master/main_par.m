function out = main_par(inputs)

% MAIN_PAR Example main program file for parallel execution.
%    Note: the name of this file is strict, owing to the constraints of 
%    parallel computing transparency.

% .. trivial example follows ..
out{1} = inputs;
out{2} = rand;
pause(rand)

