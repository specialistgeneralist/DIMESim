function output = perceived_success_after_IR_vec(u_broadcast, D_value, F_threshold)

n = length(D_value);
Z = rand(n,1);
output = ((u_broadcast == -1) .* ones(n,1) | ((100 .* ones(n,1) - D_value) .* Z > 100 .* F_threshold .* ones(n,1)));
    
output = ones(n,1) - 2.* output;
    
