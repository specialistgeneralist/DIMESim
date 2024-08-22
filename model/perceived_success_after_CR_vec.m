function output = perceived_success_after_CR_vec(S_vec, B_after_IR_succ_bool_vec, nu)

n = length(S_vec);
num_succ_neighbours_cond = (S_vec > nu .* ones(n, 1)); % boolean vector indicating if proportion of successful neighbours exceeds threshold value 
succ_IR = (B_after_IR_succ_bool_vec == 1); % boolean vector indiciating if GB reinterpreted as success after IR

output = num_succ_neighbours_cond | succ_IR;

output = ones(n,1) - 2.* output;
