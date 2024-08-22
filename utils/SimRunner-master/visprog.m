function [P] = visprog(N,NN,varargin)

% VISPROG Generate a progress bar for display.
%    P = VISPROG(N,NN) Provides a visual progress meter (string array) P given 
%    current iteration N and total iterations NN.
%
%    P = VISPROG(...,W,H) Allows for specifying the width W and height H (num 
%    lines) of the output visualisation (default W=60,H=4);
%
% See also SIMRUNNER

% Author: SA, 27 July 2005 (UNSW)

% -- default P width
if isempty(varargin)
	w = 60;
	h = 4;
else
	w = varargin{1};
	h = varargin{2};
end
dt = floor((w-4)*N/NN);
% -- create lines --
space = []; stars = [];
for i = 1:w-4
	space = [space ' '];
	stars = [stars '>'];
end
if dt > 0
	line = [' [ ' stars(1:dt) space(1:(w-4-dt)) ' ]'];
else
	line = [' [ ' space(1:w-4) ' ]'];
end
% -- construct display --
for i = 1:h
	P(i,:) = line;
end

