function [output] = perceived_success_after_IR_test(z_value, u_broadcast, D_value, F_threshold)

%Z = rand()

if (u_broadcast == -1) || ((100 - D_value) * z_value > 100 * F_threshold)
    output = -1;

else
    output = 1;

end

