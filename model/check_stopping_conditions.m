function [stopping_cond_met] = check_stopping_conditions(tt, IN)

%% Check if stopping conditions are met by assessing whether ALL mean DIME value ranges are within epsilon over a window

global D I M E A C xh x U

% Subset DIME vars to last stopping_window_steps and get mean over all agents
Ds = mean(D(:,tt - IN.stopping_window_steps + 1 - IN.rolling_average_window + 1:tt),1);
Is = mean(I(:,tt - IN.stopping_window_steps + 1 - IN.rolling_average_window + 1:tt),1);
Ms = mean(M(:,tt - IN.stopping_window_steps + 1 - IN.rolling_average_window + 1:tt),1);
Es = mean(E(:,tt - IN.stopping_window_steps + 1 - IN.rolling_average_window + 1:tt),1);

% Calculate rolling averages
D_ra = zeros(1, IN.stopping_window_steps);
I_ra = zeros(1, IN.stopping_window_steps);
M_ra = zeros(1, IN.stopping_window_steps);
E_ra = zeros(1, IN.stopping_window_steps);
for t = 1:IN.stopping_window_steps
    D_ra(t) = mean(Ds(t:t + IN.rolling_average_window - 1));
    I_ra(t) = mean(Is(t:t + IN.rolling_average_window - 1));
    M_ra(t) = mean(Ms(t:t + IN.rolling_average_window - 1));
    E_ra(t) = mean(Es(t:t + IN.rolling_average_window - 1));
end

% Calculate pc variance of DIME rolling averages as fraction of global mean for each variable over entire window
D_f = (max(D_ra) - min(D_ra)) / mean(D_ra);
I_f = (max(I_ra) - min(I_ra)) / mean(I_ra);
M_f = (max(M_ra) - min(M_ra)) / mean(M_ra);
E_f = (max(E_ra) - min(E_ra)) / mean(E_ra);
% D_f = (max(Ds) - min(Ds)) / mean(Ds);
% I_f = (max(Is) - min(Is)) / mean(Is);
% M_f = (max(Ms) - min(Ms)) / mean(Ms);
% E_f = (max(Es) - min(Es)) / mean(Es);

% Check if all DIME mean value ranges are within epsilon
if all([D_f I_f M_f E_f] < IN.stopping_epsilon)

    if IN.verbose
        fprintf('Stopping conditions met at t = %d\n', tt);
    end

    stopping_cond_met = true;

    % -- truncate all output matrices to account for stopping
    U(tt + 2 : end) = [];
    D(:, tt + 2 : end) = [];
    I(:, tt + 2: end) = [];
    M(:, tt + 2: end) = [];
    E(:, tt + 2: end) = [];
    A(:, tt + 2: end) = [];
    C(:, tt + 2: end) = [];
    xh(:, tt + 2: end) = [];
    x(:, tt + 2: end) = [];

else
    stopping_cond_met = false;

end

end

