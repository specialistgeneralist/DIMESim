% add paths
addpath(genpath(pwd))

pause on; % enables pausing
format compact; format shortG; % formatting output in command window
% dbstop if error; % enter debug mode when error arises
dbclear if error; % return to normal setting

% 'none'; 'LaTeX': interpret only math envt, i.e. $$; 'TeX': always
% interpret
set(0, 'DefaultTextInterpreter', 'TeX') %
set(0, 'DefaultLegendInterpreter', 'TeX')
set(0, 'DefaultAxesTickLabelInterpreter', 'TeX')