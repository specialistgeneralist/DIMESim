function output = action_intention(D, I, M, E)

%% If D >= (1/3) * (I + M + E), then output = 0, else 1

output = 1 - [ D >= (1/3) * (I + M + E)];

end