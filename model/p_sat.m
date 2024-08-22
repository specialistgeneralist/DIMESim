function output = p_sat(input)

% // SA: slightly different logic here to python code:
% if number_in <= 0:
%     number_out = 0
% elif number_in >= 100:
%     number_out = 100
% else:
%     number_out = number_in
% return number_out

output = input;
ix_ls0 = input <= 0;
ix_gr100 = input >= 100;
output(ix_ls0) = 0;
output(ix_gr100) = 100;

end
