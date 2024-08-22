function output = action_focus(I, M, E)

%% If I > (1/2) * (M + E) then output = -1, else output = 1

test =  I > (1/2) * (M + E);
output = zeros(size(test));
output(test)  = -1;
output(~test) = 1;

end