function InitialiseModel(IN)

%% Initialise all DIME, intention and action vars.

%% // Global vars (common) (see below for additionals based on model_variant)
global U D I M E A C xh x
global G_Outcome_D G_Outcome_I G_Outcome_M G_Outcome_E
global G_Type_D G_Type_I G_Type_M G_Type_E
global G_Interaction_D G_Interaction_I G_Interaction_M G_Interaction_E

%% // Pre-allocate
% .. prep init objects
z_1T  = zeros(1,    IN.T + 1);
z_n1  = zeros(IN.n, 1);
z_nT  = zeros(IN.n, IN.T + 1);

% -- broadcast signal
U = z_1T;      % U := global broadcast signal

% -- additional Broadcast signals/state variables as required
if ~contains(IN.model_variant, {'GeneralBroadcast'})    % -- do for everything other than General Broadcast

        global B_after_IR
        B_after_IR = z_nT; 

end
if contains(IN.model_variant, 'CollectiveReinterpretation')

    global S B_after_IR_succ_bool B_after_CR
    global adjMatrix pMatrix 

    S = z_nT;
    B_after_IR_succ_bool = z_nT;
    B_after_CR = z_nT;

    % -- set up adjMatrix
    switch IN.network
        case 'Watts-Strogatz'
            ws_graph = WattsStrogatz(IN.n, IN.k, IN.beta);
            adjMatrix = adjacency(ws_graph);
        case 'Holme-Kim'
            [~, adjMatrix] = HolmeKim(IN.n, IN.m, IN.mt, IN.N0);            
        otherwise % read adjacency matrix from file
            adjMatrix = readmatrix(IN.network);       
    end    

    % Build degree vector + matrix
    
    degVector = sum(adjMatrix, 2);

    % -- assert that the size of the matrix matches the number of agents
    assert(size(adjMatrix, 1) == IN.n, 'The size of the adjacency matrix does not match the number of agents.')

    % -- build pMatrix 
    if det(adjMatrix) == 0
        % Vectorise the below?
        for ii = 1:IN.n
            if degVector(ii) == 0
                degVector(ii) = degVector(ii) - 1;              % for agents with no neighbors, make LHS negative so that it can never exceed threshold
                adjMatrix(ii,ii) = adjMatrix(ii, ii) + 1;       % to ensure negative entry for pMatrix
            end
        end
    end
    degMatrix = diag(degVector);
    pMatrix = degMatrix^(-1) * adjMatrix;

end

% -- DIME vars
D = z_nT;     % D := Disidentification signals of all n actors
I = z_nT;     % I := Innovation signals of all n actors
M = z_nT;     % M := Moralisation signals of all n actors
E = z_nT;     % E := Energisation signals of all n actors

% -- Action and intention vars
A  = z_nT;    % A := Action intention signals of all n actors
C  = z_nT;    % C := Current action focus signals of all n actors
xh = z_nT;    % xh := Current action signals of all n actors
x  = z_nT;    % x := Action signals of all n actors

% -- Gradient vars: Outcome, Type, Interaction
G_Outcome_D = z_n1;
G_Outcome_I = z_n1;
G_Outcome_M = z_n1;
G_Outcome_E = z_n1;

G_Type_D = z_n1;
G_Type_I = z_n1;
G_Type_M = z_n1;
G_Type_E = z_n1;

G_Interaction_D = z_n1;
G_Interaction_I = z_n1;
G_Interaction_M = z_n1;
G_Interaction_E = z_n1;

%% // First time step initialise actions and action intentions --> 1

switch lower(IN.initial_action)
    case 'random' % uniformly random for each agent
        A(:, 1) = randi(2, IN.n, 1) - 1; % Action intention
        C(:, 1) = 2 * randi(2, IN.n, 1) - 3; % Current action focus
        xh(:, 1) = 2 * randi(2, IN.n, 1) - 3; % Last active action.        
    case 'all conventional'
        A(:, 1) = 1; % active
        C(:, 1) = 1; % consistent
        xh(:, 1) = 1; % conventional
    case 'all inactive conventional'
        A(:, 1) = 0; % inactive
        C(:, 1) = 1; % consistent
        xh(:, 1) = 1; % conventional
    case 'all radical'
        A(:, 1) = 1; % active
        C(:, 1) = 1; % consistent
        xh(:, 1) = -1; % radical
    case 'all inactive radical'
        A(:, 1) = 0; % inactive
        C(:, 1) = 1; % consistent
        xh(:, 1) = -1; % radical    
end
x(:, 1) = action(A(:, 1), C(:, 1), xh(:, 1));

%% // Initialise DIME vars
D(:, 1) = p_sat(normrnd(IN.studyParams.D_init_MEAN, IN.studyParams.D_init_SD, IN.n, 1));
I(:, 1) = p_sat(normrnd(IN.studyParams.I_init_MEAN, IN.studyParams.I_init_SD, IN.n, 1));
M(:, 1) = p_sat(normrnd(IN.studyParams.M_init_MEAN, IN.studyParams.M_init_SD, IN.n, 1));
E(:, 1) = p_sat(normrnd(IN.studyParams.E_init_MEAN, IN.studyParams.E_init_SD, IN.n, 1));


%% // Initialise all gradient vars.
G_Outcome_D_MEAN = IN.studyParams.DIME_grad_means(1,1);
G_Outcome_I_MEAN = IN.studyParams.DIME_grad_means(2,1);
G_Outcome_M_MEAN = IN.studyParams.DIME_grad_means(3,1);
G_Outcome_E_MEAN = IN.studyParams.DIME_grad_means(4,1);

G_Type_D_MEAN = IN.studyParams.DIME_grad_means(1,2);
G_Type_I_MEAN = IN.studyParams.DIME_grad_means(2,2);
G_Type_M_MEAN = IN.studyParams.DIME_grad_means(3,2);
G_Type_E_MEAN = IN.studyParams.DIME_grad_means(4,2);

G_Interaction_D_MEAN = IN.studyParams.DIME_grad_means(1,3);
G_Interaction_I_MEAN = IN.studyParams.DIME_grad_means(2,3);
G_Interaction_M_MEAN = IN.studyParams.DIME_grad_means(3,3);
G_Interaction_E_MEAN = IN.studyParams.DIME_grad_means(4,3);

% DIME gradient standard errors
G_Outcome_D_SE = IN.studyParams.DIME_grad_ses(1,1);
G_Outcome_I_SE = IN.studyParams.DIME_grad_ses(2,1);
G_Outcome_M_SE = IN.studyParams.DIME_grad_ses(3,1);
G_Outcome_E_SE = IN.studyParams.DIME_grad_ses(4,1);

G_Type_D_SE = IN.studyParams.DIME_grad_ses(1,2);
G_Type_I_SE = IN.studyParams.DIME_grad_ses(2,2);
G_Type_M_SE = IN.studyParams.DIME_grad_ses(3,2);
G_Type_E_SE = IN.studyParams.DIME_grad_ses(4,2);

G_Interaction_D_SE = IN.studyParams.DIME_grad_ses(1,3);
G_Interaction_I_SE = IN.studyParams.DIME_grad_ses(2,3);
G_Interaction_M_SE = IN.studyParams.DIME_grad_ses(3,3);
G_Interaction_E_SE = IN.studyParams.DIME_grad_ses(4,3);

%% // Form the DIME gradient vectors
G_Outcome_D = normrnd(G_Outcome_D_MEAN, G_Outcome_D_SE, IN.n,1);
G_Outcome_I = normrnd(G_Outcome_I_MEAN, G_Outcome_I_SE, IN.n,1);
G_Outcome_M = normrnd(G_Outcome_M_MEAN, G_Outcome_M_SE, IN.n,1);
G_Outcome_E = normrnd(G_Outcome_E_MEAN, G_Outcome_E_SE, IN.n,1);

G_Type_D = normrnd(G_Type_D_MEAN, G_Type_D_SE, IN.n,1);
G_Type_I = normrnd(G_Type_I_MEAN, G_Type_I_SE, IN.n,1);
G_Type_M = normrnd(G_Type_M_MEAN, G_Type_M_SE, IN.n,1);
G_Type_E = normrnd(G_Type_E_MEAN, G_Type_E_SE, IN.n,1);

G_Interaction_D = normrnd(G_Interaction_D_MEAN, G_Interaction_D_SE, IN.n,1);
G_Interaction_I = normrnd(G_Interaction_I_MEAN, G_Interaction_I_SE, IN.n,1);
G_Interaction_M = normrnd(G_Interaction_M_MEAN, G_Interaction_M_SE, IN.n,1);
G_Interaction_E = normrnd(G_Interaction_E_MEAN, G_Interaction_E_SE, IN.n,1);
