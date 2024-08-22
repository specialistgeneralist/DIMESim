function [x_con, x_rad, x_ina, xi_conhab, xi_coninv, xi_radhab, xi_radinv, xi_conlat, xi_radlat, x_conhab, x_coninv, x_radhab, x_radinv] = std_outputs_actions(A, C, xh, x)
%{
x_con, x_rad, x_ina = Timeseries of Conventional, Radical, Inactive fractions
xi_conhab, xi_coninv, = Timeseries of inactive fractions whose xh is conventional and
current focus is conventional (conhab), radical (coninv)
xi_conlat, xi_radlat = Timeseries of inactive fractions whose
current focus is conventional (conlat), radical (radlat)
%}

epsilon = 0;

% Calculate the number of actors and the number of time-steps in the simulation
num_actors_scalar = size(A, 1);
num_time_steps = size(A,2); %size(time_set, 1); time_set is currently row vector. changing to number of columns in A to avoid ambiguity in time_set shapes, and since later matrix calculations require this value anyway

% Calculate the number of active and inactive actors at each time-step
num_active_actors = ones(1, num_actors_scalar) * A;
num_actors_vector = num_actors_scalar .* ones(1, num_time_steps);
num_inactive_actors = num_actors_vector - num_active_actors;

% Calculate total numbers of conventional, radical actors at each time-step
% Radical actors
ones_array = ones(size(x));
temp = (1/2) .* (ones_array - x);                   % eliminating conventional actors from count
temp = temp + epsilon .* ones_array;
temp = floor(temp);                                 % eliminating inactive actors from count
num_rad_actors = ones(1, num_actors_scalar) * temp;

% Conventional actors
num_con_actors = num_active_actors - num_rad_actors;

% Conventional, radical and inactive percentages
x_con = (1/num_actors_scalar) .* num_con_actors;
x_rad = (1/num_actors_scalar) .* num_rad_actors;
x_ina = (1/num_actors_scalar) .* num_inactive_actors;

% Habitual and innovative active actors
x_conhab = sum(A == 1 & C == 1 & xh == 1) / num_actors_scalar;
x_coninv = sum(A == 1 & C == -1 & xh == -1) / num_actors_scalar;
x_radhab = sum(A == 1 & C == 1 & xh == -1) / num_actors_scalar;
x_radinv = sum(A == 1 & C == -1 & xh == 1) / num_actors_scalar;

% Classify inactive actor states

% 1) Determine inactive actors with "conventional" and "radical" xh values 
inactive_actors = ones(size(A)) - A;   % Array of all the inactive actors over the simulation

% Get the historical action of only the inactive actors (Hadamard product
% of xh and inactive actors array)
xh_inactive = xh .* inactive_actors;

% Extract the number of inactive actors who display radical xh at each
% time-step
temp = (1/2) .* (ones_array - xh_inactive);         % eliminating inactive conventional xh from count
temp = temp + epsilon .* ones_array;
temp = floor(temp);                                 % eliminating active actors from count
inactive_radxh_actors = temp;
num_inactive_radxh_actors = ones(1, num_actors_scalar) * inactive_radxh_actors;

% Since xh can only be either "conventional" or "radical", the number of
% inactives with conventional xh is thus determined by:
num_inactive_conxh_actors = num_inactive_actors - num_inactive_radxh_actors;

% The array which identifies the actors with conventional xh is thus
inactive_conxh_actors = inactive_actors - inactive_radxh_actors;

% 2) Determine the (xh,C) state of each inactive actors by taking the
% Hadamard products of the xh arrays (produced above) with the C array

inactive_radxh_curr_focus = inactive_radxh_actors .* C;

% The entries of this matrix which are '1' denote actors who are 'radical 
% habitual' (radhab); the entries which are '-1' denote actors who are '
% radical innovative' (radinv)

% Extract the radinv actors (-1 entries) using the mathod as before
temp = (1/2) .* (ones_array - inactive_radxh_curr_focus);
temp = temp + epsilon .* ones_array;
temp = floor(temp);
num_xi_radinv = ones(1, num_actors_scalar) * temp;

% Since each actors can only have either C = 1 or C = -1, we find the
% number of inactive actors with radical xh and C = 1 (habitual) by

num_xi_radhab = num_inactive_radxh_actors - num_xi_radinv;

% We do the same for the inactive actors with conventional xh
inactive_conxh_curr_focus = inactive_conxh_actors .* C;

temp = (1/2) .* (ones_array - inactive_conxh_curr_focus);
temp = temp + epsilon .* ones_array;
temp = floor(temp);
num_xi_coninv = ones(1, num_actors_scalar) * temp;

num_xi_conhab = num_inactive_conxh_actors - num_xi_coninv;


% Calculate proportions
xi_conhab = (1/num_actors_scalar) * num_xi_conhab;
xi_coninv = (1/num_actors_scalar) * num_xi_coninv;
xi_radhab = (1/num_actors_scalar) * num_xi_radhab;
xi_radinv = (1/num_actors_scalar) * num_xi_radinv;

% Latent actors are those inactive who possess this focus currently
xi_conlat = xi_conhab + xi_radinv;
xi_radlat = xi_coninv + xi_radhab;

end