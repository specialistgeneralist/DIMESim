function res_contourf(RESFNAME, PAR_VARX, PAR_VARY, varargin)

%% Analyse final outcomes after varying two parameters (x and y)

% -- constant settings
in = struct('RESFNAME', RESFNAME);
in.ps = plot_settings_common('contourf');
dime_set_full = {'D_avg' 'I_avg' 'M_avg' 'E_avg'};
% action_set_full = {'Con_f', 'Rad_f', 'InaConLat_f', 'InaRadLat_f', 'DominantAction', 'DominantAction3D'};
action_set_full = {'ConHab_f', 'Inv_f', 'RadHab_f', 'InaConHab_f', 'InaInv_f', 'InaRadHab_f', 'DominantAction', 'DominantAction3D'};

% -- input parameters
PAPER_FIGURE = []; % -- (3--5, -1 -- -4) Number of relevant figure from paper. Loads hardcoded settings for all parameters

COMPARTMENTALIZE = false;   % -- whether to move input and output files related to this run into a separate folder
in.FIGURE_SIZE = [];% -- figure window size
INCLUDE_LABELS = true;           % -- whether to include axis, tick, contour labels
in.LIGHTNESS_MAX = 1;      % maximum lightness for color range, in fraction. 1 = white.
% in.MARK_UNCONVERGED = false;     % -- whether to mark unconverged points
in.MINIMAL_DETAIL = false;  % -- false for sufficient detail for standalone
in.N_COLOURS = 256;             % -- number of colors to use in heatmap plot. Use 2 for low/high contrast
PAR_CONST_SET = [];     % -- string ("") array of names of constant parameters to be printed in figure. e.g. ["model_variant", "nu", "R"]
PAR_SUBSET = struct();  % -- struct where field is parameter and value is list of parameter values used to filter data set
PLOT_GROUP_BY = 'Metric'; % -- type of collection plots will be gathered into. Individual - every plot is on a separate figure. Parameters - subplots share same set of parameters; Metric - subplots share same metric being plotted
PLOT_RC = []; % number of rows and columns for subplot arrangement
PRINT_PLOT = '';        % -- print format for plots. if '', then don't print
STEM_SUBSET = {};       % -- cell array subset of DIME and action outputs to be plotted. Default: empty, print all
SUBPLOT_SUBTITLE = false;   % -- whether titles of subplots should contain subtitle added to full figure (for cropping to subplot)

% -- handle varargin as name-value pairs. Name is not case-sensitive
for i = 1:2:length(varargin)
    NAME = lower(varargin{i});
    VALUE = varargin{i+1};
    switch NAME
        case 'compartmentalize'
            COMPARTMENTALIZE = VALUE;
        case 'figure_size'
            in.FIGURE_SIZE = VALUE;
        case 'include_labels'
            INCLUDE_LABELS = VALUE;
        case 'lightness_max'
            in.LIGHTNESS_MAX = VALUE;
        % case 'mark_unconverged'
        %     in.MARK_UNCONVERGED = VALUE;
        case 'minimal_detail'
            in.MINIMAL_DETAIL = VALUE;
        case 'n_colours'
            in.N_COLOURS = VALUE;
        case 'paper_figure'
            PAPER_FIGURE = VALUE;
        case 'par_const_set'
            PAR_CONST_SET = VALUE;
        case 'par_subset'
            PAR_SUBSET = VALUE;
        case 'plot_group_by'
            PLOT_GROUP_BY = VALUE;
        case 'plot_rc'
            PLOT_RC = VALUE;
        case 'print_plot'
            PRINT_PLOT = VALUE;
        case 'stem_subset'
            STEM_SUBSET = VALUE;
        case 'subplot_subtitle'
            SUBPLOT_SUBTITLE = VALUE;            
    end
end

% -- paper figure settings
switch PAPER_FIGURE
    case 3
        RESFNAME = 'ParameterSweep_p-F.mat';
        PAR_VARX = 'p';
        PAR_VARY = 'F';
        PLOT_GROUP_BY = 'Parameters';
        STEM_SUBSET = {'ConHab_f', 'Inv_f', 'RadHab_f', 'InaConHab_f', 'InaInv_f', 'InaRadHab_f'};
        in.MINIMAL_DETAIL = true;
        in.FIGURE_SIZE = [1097 889];    
    case 4
        RESFNAME = 'ParameterSweep_p-F.mat';
        PAR_VARX = 'p';
        PAR_VARY = 'F';
        PLOT_GROUP_BY = 'DominantAction';       
        in.MINIMAL_DETAIL = true;
        in.FIGURE_SIZE = [700 465];
    case 5
        RESFNAME = 'ParameterSweep_nu-R_p-F.mat';
        PAR_VARX = 'p';
        PAR_VARY = 'F';
        PAR_SUBSET = struct('nu', [], 'R', []);
        PLOT_GROUP_BY = 'Metric';
        STEM_SUBSET = {'DominantAction3D'};
        in.MINIMAL_DETAIL = true;
        in.FIGURE_SIZE = [];%[1097 889];
    % Supplementary Section
    case -1
        RESFNAME = 'ParameterSweep_p-F.mat';
        PAR_VARX = 'p';
        PAR_VARY = 'F';
        PLOT_GROUP_BY = 'Parameters';
        STEM_SUBSET = {'D_avg' 'I_avg' 'M_avg' 'E_avg'};
        in.MINIMAL_DETAIL = true;
        in.FIGURE_SIZE = [1097 369];
    case -2
        RESFNAME = 'ParameterSweep_nu-R_p-F.mat';
        PAR_VARX = 'p';
        PAR_VARY = 'F';
        PAR_SUBSET = struct('nu', [], 'R', []);
        PLOT_GROUP_BY = 'Metric';
        STEM_SUBSET = {'ConHab_f', 'InaConHab_f', 'Inv_f', 'InaInv_f', 'RadHab_f', 'InaRadHab_f', 'DominantAction', 'D_avg', 'I_avg', 'M_avg', 'E_avg'};        
        in.MINIMAL_DETAIL = true;
        in.FIGURE_SIZE = [];%[1097 889];
    case -3
        RESFNAME = 'ParameterSweep_nu-R.mat';
        PAR_VARX = 'nu';
        PAR_VARY = 'R';
        PLOT_GROUP_BY = 'Parameters';
        STEM_SUBSET = {'ConHab_f', 'Inv_f', 'RadHab_f', 'InaConHab_f', 'InaInv_f', 'InaRadHab_f'};
        in.MINIMAL_DETAIL = true;
        in.FIGURE_SIZE = [1097 889];
    case -4
        RESFNAME = 'ParameterSweep_nu-R.mat';
        PAR_VARX = 'nu';
        PAR_VARY = 'R';
        PLOT_GROUP_BY = 'Parameters';
        STEM_SUBSET = {'D_avg', 'I_avg', 'M_avg', 'E_avg'};
        in.MINIMAL_DETAIL = true;
        in.FIGURE_SIZE = [1097 889];
end

% -- load results
RESFNAME_var = who('-file', RESFNAME);
RESFNAME_split = split(RESFNAME, {filesep, '.mat'});
simulation_name = RESFNAME_split{end-1};
if ismember('summary', RESFNAME_var)
    load(RESFNAME, 'summary');
    ALPHA = 0.05;           % -- confidence interval alpha
    NBOOT = 1000;           % -- number of bootstrap samples to use for confidence intervals
    for i = 1:numel(summary)
        summary(i).inputs.alpha = ALPHA;
        summary(i).inputs.nboot = NBOOT;
    end

    [M, S, inputs_collated] = average_results(summary, {PAR_VAR});
else
    load(RESFNAME, 'M', 'S', 'inputs_collated');
end

% -- filtering data set
PAR_SUBSET_fields = fields(PAR_SUBSET);
n_ps = numel(PAR_SUBSET_fields);
PAR_SUBSET_array = repmat(struct('Name', '', 'Values', []), n_ps, 1); % restructing as array for easier referencing
Mrow_keep = true(size(M,1),1); % rows to keep
for k =1:n_ps
    par = PAR_SUBSET_fields{k};
    par_val = PAR_SUBSET.(par);
    if ~isempty(par_val)
        Mrow_keep = Mrow_keep & ismember(M.(par), par_val);
    else % if empty, use all values
        par_val = unique(M.(par))';
        PAR_SUBSET.(par) = par_val;
    end
    PAR_SUBSET_array(k) = struct('Name', par, 'Values', par_val);
    inputs_collated.(par) = par_val;
end
M = M(Mrow_keep, :);
S = S(Mrow_keep, :);

    %% Plot Settings
% -- creating plot subtitle with parameter settings
group_subtitle = ""; % enclosing in $$ for LaTeX interpretation
if INCLUDE_LABELS    
    for par = PAR_CONST_SET
        par_value = inputs_collated.(par);
        switch par
            case 'model_variant'
                par_string = string(par_value(isstrprop(par_value,'upper')));     % using abbreviation for model variant
            case 'nu'
                par_string = append("\phi=", string(par_value));
            case 'n_replicates'
                par_string = append("n_{replicates}=", string(par_value));
            otherwise
                par_string = append(string(par), "=",string(par_value));
        end
        if numel(par_string) > 1
            par_string = append("[", strjoin(par_string, ','), "]");
        end

        if ~isequal(PLOT_GROUP_BY, 'Metric') || ~ismember(par, PAR_SUBSET_fields)
            group_subtitle = append(group_subtitle, par_string, '; ');        
        end
    end
end
if isempty(STEM_SUBSET)
    dime_set = dime_set_full;
    action_set = action_set_full;    
else
    dime_set = dime_set_full(ismember(dime_set_full, STEM_SUBSET));
    action_set = action_set_full(ismember(action_set_full, STEM_SUBSET));    
end

% -- Plot :: heatmap
filename_prefix = erase(RESFNAME, '.mat');
if INCLUDE_LABELS && (SUBPLOT_SUBTITLE || ~ismember(PLOT_GROUP_BY, {'Parameter', 'Metric'}))
    in.plot_subtitle = group_subtitle;
else
    in.plot_subtitle = '';
end
% in.action_stems = {'Con_f', 'Rad_f', 'InaConLat_f', 'InaRadLat_f'};
in.action_stems = {'ConHab_f', 'RadHab_f', 'Inv_f', 'InaConHab_f', 'InaRadHab_f', 'InaInv_f'}; % for determining dominant action 
in.group_subtitle = group_subtitle;
in.PAR_VARX = PAR_VARX;
in.PAR_VARY = PAR_VARY;
in.xticks = unique(M.(PAR_VARX));
in.yticks = unique(M.(PAR_VARY));
in.nx = numel(in.xticks); % number of elements PAR_VARX was varied over
in.ny = numel(in.yticks); % number of elements PAR_VARY was varied over
in.INCLUDE_LABELS = INCLUDE_LABELS;
in.PLOT_GROUP_BY = PLOT_GROUP_BY;
in.PAR_SUBSET_array = PAR_SUBSET_array;
in.n_ps = n_ps;

    %% Plotting
switch PLOT_GROUP_BY    
    case 'Individual'
        fi = 0;
        in.compact_axis_label = false;

        plot_zlim = [0 100];
        zscaling = 1;
        for k = 1:numel(dime_set)
            fi = fi + 1;            
            this_var = dime_set{k};            
            fh = plot_single(fi, this_var, M, plot_zlim, zscaling, in);
            stem_ID = extractBefore([this_var '_'],'_');
            fh.Name = [simulation_name '_' stem_ID];            
        end
        
        plot_zlim = [0 100];
        zscaling = 100;
        for k = 1:numel(action_set)
            fi = fi + 1;            
            this_var = action_set{k};
            fh = plot_single(fi, this_var, M, plot_zlim, zscaling, in);            
            stem_ID = extractBefore([this_var '_'],'_');
            fh.Name = [simulation_name '_' stem_ID];
        end
       
    case 'Parameters'
        in.compact_axis_label = false;
        
        if ~isempty(dime_set)
            output_category = "DIME Value"; 
            plot_zlim = [0 100];
            zscaling = 1;
            fh = plot_by_inputgroup(1, dime_set, M, [1 numel(dime_set)], plot_zlim, zscaling, output_category, in);
            fh.Name = [simulation_name '_DIME'];
            if string(PRINT_PLOT) == "png"
                print(fh, '-dpng', '-r100', [filename_prefix '_dime.png'])
            end
        end
        
        if ~isempty(action_set)
            output_category = "Action Percentage"; 
            plot_zlim = [0 100];
            zscaling = 100;
            fh = plot_by_inputgroup(2, action_set, M, [2 ceil(numel(action_set)/2)], plot_zlim, zscaling, output_category, in);
            fh.Name = [simulation_name '_Action'];
            if string(PRINT_PLOT) == "png"
                print(fh, '-dpng', '-r100', [filename_prefix '_actions.png'])
            end
        end

    case 'Metric'
        fi = 0;
        in.compact_axis_label = true;
        in.plot_rc = [1 1];        
        for k = 1:n_ps
            in.plot_rc(3 - k) = numel(PAR_SUBSET_array(k).Values); % reverse order, columns then rows        
        end        
        
        plot_zlim = [0 100];
        zscaling = 1;
        for k = 1:numel(dime_set)
            fi = fi + 1;
            this_var = dime_set{k};
            stem_ID = extractBefore([this_var '_'],'_');            
            fh = plot_by_outputgroup(fi, this_var, M, plot_zlim, zscaling, in);
            fh.Name = [simulation_name '_' stem_ID];
        end

        plot_zlim = [0 100];
        zscaling = 100;
        for k = 1:numel(action_set)
            fi = fi + 1;
            this_var = action_set{k};
            stem_ID = extractBefore([this_var '_'],'_');
            fh = plot_by_outputgroup(fi, this_var, M, plot_zlim, zscaling, in);
            fh.Name = [simulation_name '_' stem_ID];
        end

    case 'DominantAction'
        in.compact_axis_label = false;
        output_category = "DominantAction Percentage"; 
        plot_zlim = [0 100];
        zscaling = 100;
        fh = plot_by_inputgroup(1, {'DominantAction', 'DominantAction3D'}, M, [1 2], plot_zlim, zscaling, output_category, in);
        fh.Name = [simulation_name '_DominantAction'];
        if string(PRINT_PLOT) == "png"
            print(fh, '-dpng', '-r100', [filename_prefix '_actions.png'])
        end        

    case 'Average Action'        
        in.compact_axis_label = false;
        in.stems = {'Con_f', 'Rad_f', 'InaConLat_f', 'InaRadLat_f'};        
        plot_zlim = [0 100];
        zscaling = 100;
               
        fh = figure_subplot(1); clf;
        fh.Name = [simulation_name '_AverageAction'];
        plot_AverageAction(fh, M, plot_zlim, zscaling, in);      
end

    %% Post-Plot
if COMPARTMENTALIZE
    fold = filename_prefix;
    mkdir(fold);
    addpath(fold);
    copyfile([filename_prefix '.txt'], fold);
    movefile(RESFNAME, fold);
    movefile([filename_prefix '_dime.png'], fold);
    movefile([filename_prefix '_actions.png'], fold);
end

end



% --------------------------
function fh = plot_by_inputgroup(fi, var_stems, M, plot_rc, zlim, zscaling, output_category, in)
%%% plots subplots that have same input (parameters) in same figure

fh = figure_subplot(fi, plot_rc, in.FIGURE_SIZE); % 1.01*[465 520] .* flip(plot_rc)
tl = tiledlayout(plot_rc(1), plot_rc(2), 'Padding', 'Compact');
for k = 1:numel(var_stems)    
    h = nexttile;
    this_var = var_stems{k};        
    plot_single(h, this_var, M, zlim, zscaling, in);
end
if in.INCLUDE_LABELS
    plot_sgtitle = output_category + " (end of exp): " + in.group_subtitle;
else
    plot_sgtitle = output_category;
end
% labelwrap(fh, plot_sgtitle, 'sgtitle');
if ~in.MINIMAL_DETAIL
    title(tl, plot_sgtitle,'Interpreter', 'tex'); % labelwrap has issues. perhaps due to asking for dimensions of tiledlayout
end
end



% --------------------------
function fh = plot_by_outputgroup(fi, this_var, M, zlim, zscaling, in)
%%% plots subplots that have same output (action or dime) in same figure

stem_ID = extractBefore([this_var '_'],'_');
psa = in.PAR_SUBSET_array;

fh = figure_subplot(fi, in.plot_rc, [465 520] .* in.plot_rc);
tl = tiledlayout(in.plot_rc(1), in.plot_rc(2), 'Padding', 'Tight');
if in.n_ps >= 1    
    par1 = psa(1).Name;
    if par1 == "nu"
        par1_subtitlename = "\phi";
    else
        par1_subtitlename = par1;
    end
    if in.n_ps == 2    
        par2 = psa(2).Name;
        if par2 == "nu"
            par2_subtitlename = "\phi";
        else
            par2_subtitlename = par2;
        end
    end
end

for r = in.plot_rc(1):-1:1 % plotting starts from toprow
    for c = 1:in.plot_rc(2)        
        h = nexttile;
                
        if in.n_ps >= 1
            par1_Val = psa(1).Values(c);
            if in.n_ps < 2
                Mrow_plot = M.(par1) == par1_Val;
                in.plot_subtitle = append(par1_subtitlename, '=', string(par1_Val));
            else
                par2_Val = psa(2).Values(r);
                Mrow_plot = M.(par1) == par1_Val & M.(par2) == par2_Val;
                in.plot_subtitle = append(par1_subtitlename, '=', string(par1_Val) ...
                    , '; ', par2_subtitlename, '=', string(par2_Val));
            end
        else
            Mrow_plot = true(size(M,1),1); % rows to keep
        end            
        M_plot = M(Mrow_plot, :);        
        [~, cb] = plot_single(h, this_var, M_plot, zlim, zscaling, in);        

        if h.Layout.Tile == 1 %r == 1 && c == in.plot_rc(2) % use colorbar of first tile for full layout            
            if ismember(this_var, {'DominantAction', 'DominantAction3D'})
                cb.Layout.Tile = 'North';
            else
                cb.Layout.Tile = 'East';
            end
        else
            colorbar('off'); % turn off other colorbars
        end
    end
end
plot_sgtitle = [in.ps.(stem_ID).DisplayName ' (end of exp)'];
if ~in.MINIMAL_DETAIL    
    title(tl, {plot_sgtitle, in.group_subtitle}); % labelwrap has issues. perhaps due to asking for dimensions of tiledlayout
end
if in.n_ps >= 1        
    xlabel(tl, in.ps.(par1).DisplayName);        
    if in.n_ps == 2            
        ylabel(tl, in.ps.(par2).DisplayName);
    end
end


end



% --------------------------
function [fh, cb] = plot_single(fi, this_var, M, plot_zlim, zscaling, in)
    switch this_var
        case 'DominantAction'            
            fh = figure_subplot(fi, [], [465 520]);
            cb = plot_DominantAction3D(fh, M, plot_zlim, zscaling, in);
            view(2);
        case 'DominantAction3D'
            fh = figure_subplot(fi, [], [650 780]);            
            cb = plot_DominantAction3D(fh, M, plot_zlim, zscaling, in);
            view(-45.38, 33.184); % custom fine-tuning
            daspect([1 1 90]);
            align_axislabel([], gca);
        otherwise
            fh = figure_subplot(fi, [], [465 520]);
            cb = plot_regular(fh, this_var, M, plot_zlim, zscaling, in);
    end
end


% --------------------------
function cb = plot_regular(h, this_var, M, plot_zlim, zscaling, in)
%%% plots single contourf

hold on;
PAR_VARX = in.PAR_VARX;
PAR_VARY = in.PAR_VARY;
stem_ID = extractBefore(this_var,'_');

z_values = M.(this_var) * zscaling;
z_range = round(linspace(plot_zlim(1), plot_zlim(2), 11)); % for contour levels
while sum(z_values < z_range(2)) > numel(z_values)/3 % if at least one third of values are below second-lowest contour
    % z_range = [z_range(1) round(mean(z_range(1:2)),2,'significant') z_range(2:end)]; % insert half-way contour. rounding new value to 2 significant digits
    z_range = [z_range(1) mean(z_range(1:2)) z_range(2:end)]; % insert half-way contour. 
end
z_range = round(z_range,2,'significant'); % rounding all countours to 2 significant digits
[~, ch] = contourf(reshape(M.(PAR_VARX),in.nx,in.ny), ...
             reshape(M.(PAR_VARY),in.nx,in.ny), ...
             reshape(z_values, in.nx,in.ny), ...
             z_range, "ShowText",in.INCLUDE_LABELS); %"FaceAlpha",0.25
ch.LineStyle = in.ps.(stem_ID).LineStyle;
%heatmap(M, in.PAR_VARX, in.PAR_VARY, 'ColorVariable', this_var, 'ColorMap', sky(in.N_COLOURS));    
c_max = in.ps.(stem_ID).Color;
c_map = 1 - (1 - in.LIGHTNESS_MAX * (linspace(1,0,in.N_COLOURS))') * (1 - c_max); % range of base color with lightness varying from base color to max lightness parameter
colormap(h, c_map);    
set(gca, 'CLim', plot_zlim); 
axis square;
% set(gca, 'FontSize',12);

% xL = xlim();
% yL = ylim();
% xt = linspace(xL(1),xL(2),11);
% yt = linspace(yL(1),yL(2),11);
xt = in.xticks;
yt = in.yticks;
if in.ps.(stem_ID).MarkerSize > 0       % add dots overlay
    x_dots = mean([xt(1:end-1) ; xt(2:end)]);
    y_dots = mean([yt(1:end-1) ; yt(2:end)]);
    [x_dots_mesh, y_dots_mesh] = meshgrid(x_dots, y_dots);
    plot(x_dots_mesh(:), y_dots_mesh(:), 'k.', 'MarkerSize', in.ps.(stem_ID).MarkerSize);
end

plot_title = '';
switch in.PLOT_GROUP_BY
    case {'Metric', 'DominantAction'}
        plot_title = [plot_title char('a' - 1 +tilenum(h)) ') '];
    case {'Individual', 'Parameters'}
        plot_title = [plot_title in.ps.(stem_ID).DisplayName];
end
plot_title = [plot_title in.plot_subtitle];
% labelwrap(h, plot_title, 'title', 'tex');

if in.INCLUDE_LABELS
    if in.compact_axis_label
        xlabel(PAR_VARX);
        ylabel(PAR_VARY);
    else
        labelwrap(h, in.ps.(PAR_VARX).DisplayName, 'xlabel', 'tex');
        labelwrap(h, in.ps.(PAR_VARY).DisplayName, 'ylabel', 'tex');
    end
    xticks(xt);
    yticks(yt);
    cb = colorbar('Location', 'northoutside');

    switch in.PLOT_GROUP_BY
        case {'Metric', 'DominantAction'}
            title(plot_title);
        case {'Individual', 'Parameters'}
            cb.Label.String = plot_title;
    end    
    grid on
else                
    xticks([])
    yticks([])
    cb = colorbar('Ticks', plot_zlim(1) + [0.25 0.75]*(plot_zlim(2) - plot_zlim(1)), ...
        'TickLabels', ["Low" "High"]);    
    % colorbar('Ticks', []);
end
hold off;
end




% --------------------------
function cb = plot_DominantAction3D(h, M, plot_zlim, zscaling, in)
%%% plots custom heatmap with each cell colored according to dominant action

hold on;
PAR_VARX = in.PAR_VARX;
PAR_VARY = in.PAR_VARY;
c_map = [in.ps.Con.Color ; in.ps.Inv.Color; in.ps.Rad.Color];
colormap(h, c_map);

xt = unique(M.(PAR_VARX));
yt = unique(M.(PAR_VARY));
dx = xt(2) - xt(1);
dy = yt(2) - yt(1);
[xMesh, yMesh] = meshgrid(xt, yt);
zMesh = nan(size(xMesh));
cMesh = nan(size(xMesh));
actionMesh = strings(size(xMesh));

M_actions = table2array(M(:, in.action_stems));
M.MaxAction = strings(height(M),1); % type of action that is max
M.MaxActionValue = zeros(height(M),1); % value of max action
M.MaxActionColorValue = zeros(height(M),1); % color value for max action according to cmap
for i = 1:height(M_actions) %    
    action_max = max(M_actions(i,:));
    M.MaxActionValue(i) = action_max;
    stem = in.action_stems{find(M_actions(i,:) == action_max, 1)}; % find dominant action
    M.MaxAction(i) = stem;
    % M.MaxActionColorValue(i) = find(ismember(in.stems, stem)) - 1 + action_max; % color intensity varies with magnitude
    if contains(stem, 'Con')
        M.MaxActionColorValue(i) = 1;
    elseif contains(stem, 'Inv')
        M.MaxActionColorValue(i) = 2;
    elseif contains(stem, 'Rad')
        M.MaxActionColorValue(i) = 3;
    else
        M.MaxActionColorValue(i) = NaN;
        disp('Color Error');
    end

    x = find(xt == M.(PAR_VARX)(i));
    y = find(yt == M.(PAR_VARY)(i));
    zMesh(x, y) = action_max * zscaling;
    cMesh(x, y) = M.MaxActionColorValue(i);
    actionMesh(x, y) = stem;
end

b = bar3(zMesh', 1); % use 100% of allocated space; i.e. no gap between bars. bar3 plots (y,x)
view(3);
for x = 1:numel(xt)
    for y =1:numel(yt)
        rows = (y-1)*6 + (1:6);        
        b(x).XData(rows, :) = (b(x).XData(rows, :) - x) / 1 * dx + xt(x); % preserve same relative scaling
        b(x).YData(rows, :) = (b(x).YData(rows, :) - y) / 1 * dy + yt(y); % preserve same relative scaling
        b(x).CData(rows, :) = cMesh(x,y);  

        ps_stem = in.ps.(extractBefore(actionMesh(x,y),'_'));
        sh = surface(b(x).XData(rows, :), b(x).YData(rows, :), b(x).ZData(rows, :) ...
            , 'facecolor', c_map(cMesh(x,y),:), 'LineStyle', ps_stem.HatchLineStyle, 'LineWidth', ps_stem.LineWidth);
        % hatchfill2(sh, ps_stem.HatchStyle, 'HatchAngle', ps_stem.HatchAngle...
        %     , 'HatchColor', ps_stem.HatchColor, 'HatchLineStyle', ps_stem.HatchLineStyle);   
        if strcmpi(ps_stem.HatchStyle, 'single')
            plot3([xt(x) - dx/2 ; xt(x) + dx/2], [yt(y) - dy/2 ; yt(y) + dy/2], [zMesh(x,y) ; zMesh(x,y)], '--', 'Color', ps_stem.HatchColor);
        end
    end
    delete(b(x));
end

set(gca, 'CLim', [1 3], 'FontSize',10);
axis square;
% xlim([min(xt) max(xt)]);
% ylim([min(yt) max(yt)]);
zlim(plot_zlim);
% set(gca, 'FontSize',12, 'ColorLimits',in.ylim, 'Interpreter', 'None');
grid on;

plot_title = '';
if ismember(in.PLOT_GROUP_BY, {'Parameters', 'Metric', 'DominantAction'})
    plot_title = [plot_title char('a' - 1 +tilenum(h)) ') '];
end
if ismember(in.PLOT_GROUP_BY, {'Individual', 'Parameters'})
    plot_title = [plot_title in.ps.DominantAction.DisplayName];    
end
plot_title = [plot_title in.plot_subtitle];

if in.INCLUDE_LABELS
    if in.compact_axis_label
        xlabel(PAR_VARX);
        ylabel(PAR_VARY);
        zlabel('Population Fraction')
    else
        % xlabel(in.ps.(PAR_VARX).DisplayName);
        % ylabel(in.ps.(PAR_VARY).DisplayName);
        zlabel('Population Fraction')
        labelwrap(h, in.ps.(PAR_VARX).DisplayName, 'xlabel', 'tex');
        labelwrap(h, in.ps.(PAR_VARY).DisplayName, 'ylabel', 'tex');
        % labelwrap(h, in.ps.MaxAction.DisplayName, 'zlabel', 'tex');
    end
    xticks(xt);
    yticks(yt);    
    cb = colorbar('Ticks', 1 + ((1:3)-1)*2/3 + 2/6, 'TickLabels', ["Conventional" "Innovative" "Radical"], 'Location', 'northoutside');    
else                
    xticks([])
    yticks([])
    cb = colorbar('Ticks', []);
end

if strcmpi(in.PLOT_GROUP_BY, 'Metric')
    title(plot_title, 'interpreter', 'tex');
else
    cb.Label.String = plot_title; %in.ps.DominantAction.DisplayName;
end

hold off;
end



% --------------------------
function plot_AverageAction(h, M, zlim, zscaling, in)
%%% plots custom heatmap with each cell colored according to dominant action

hold on;
PAR_VARX = in.PAR_VARX;
PAR_VARY = in.PAR_VARY;
c_func = @(lightness, c_max) 1 - (1 -  lightness) * (1 - c_max); % 0 <=z<=1, z should be column vector indicating lightness of colour, c_max is row RGB vector
lightness_range = linspace(in.LIGHTNESS_MAX, 0, in.N_COLOURS);
c_map = [c_func(flip(lightness_range(2:end))', in.ps.Rad.Color) ...
    ; [1 1 1] ...
    ; c_func(lightness_range(2:end)', in.ps.Con.Color)]; % starting from 2 to avoid double-counting white
colormap(h, c_map);

xt = unique(M.(PAR_VARX));
yt = unique(M.(PAR_VARY));
dx = xt(2) - xt(1);
dy = yt(2) - yt(1);
[xMesh, yMesh] = meshgrid(xt, yt);
xPatch = [xMesh(:)' - dx/2 ; xMesh(:)' + dx/2 ; xMesh(:)' + dx/2 ; xMesh(:)' - dx/2]; % column is a point, row is a vertix of square centred at point, going anticlockwise from SW corner
yPatch = [yMesh(:)' - dy/2 ; yMesh(:)' - dy/2 ; yMesh(:)' + dy/2 ; yMesh(:)' + dy/2];
M.AverageOrientation = M.Con_f + M.InaConLat_f - M.Rad_f - M.InaRadLat_f; % -1 (100% radical) to 0 (50% rad, 50% con) to 1 (100% con)

patch(xPatch, yPatch, M.AverageOrientation, 'EdgeColor', 'black'); %, 'LineWidth', 1);
%1 - (1 - in.LIGHTNESS_MAX * (linspace(1,0,in.N_COLOURS))') * (1 - c_max); % range of base color with lightness varying from base color to max lightness parameter
M_actions = table2array(M(:, in.stems));
for i = 1:height(M_actions) % add dots overlay   
    stem = in.stems{find(M_actions(i,:) == max(M_actions(i,:)), 1)}; % find dominant action
    stem_ID = extractBefore(stem,'_');
    ms = in.ps.(stem_ID).MarkerSize;
    if ms > 0
        plot(M.(PAR_VARX)(i), M.(PAR_VARY)(i), 'k.', 'MarkerSize', ms);
    end
end

set(gca, 'CLim', [-1 1]); 
axis square;
xlim([min(xPatch(:)) max(xPatch(:))]);
ylim([min(yPatch(:)) max(yPatch(:))]);
% set(gca, 'FontSize',12, 'ColorLimits',in.ylim, 'Interpreter', 'None');

labelwrap(h, 'Average Action (end of exp)', 'title');
labelwrap(h, in.plot_subtitle, 'subtitle', 'tex');
if in.INCLUDE_LABELS
    if in.compact_axis_label
        xlabel(PAR_VARX);
        ylabel(PAR_VARY);
    else
        labelwrap(h, in.ps.(PAR_VARX).DisplayName, 'xlabel', 'tex');
        labelwrap(h, in.ps.(PAR_VARY).DisplayName, 'ylabel', 'tex');
    end
    xticks(xt);
    yticks(yt);
    colorbar('Ticks', [-0.5 0.5], 'TickLabels', ["Radical" "Conventional"]);    
else                
    xticks([])
    yticks([])
   colorbar('Ticks', []);
end
hold off;
end




% --------------------------
function cb = plot_DominantAction(h, M, zlim, zscaling, in)
%%% plots custom heatmap with each cell colored according to dominant action

hold on;
PAR_VARX = in.PAR_VARX;
PAR_VARY = in.PAR_VARY;
%{
c_func = @(lightness, c_max) 1 - (1 -  lightness) * (1 - c_max); % 0 <=z<=1, z should be column vector indicating lightness of colour, c_max is row RGB vector
lightness_range = linspace(in.LIGHTNESS_MAX, 0, in.N_COLOURS);
c_map = [];
for stem = in.stems
    stem_ID = extractBefore(stem{1},'_');
    c_map = [c_map ; c_func(lightness_range', in.ps.(stem_ID).Color)]; % compiling ranges of all colours
end
%}
c_map = [in.ps.Rad.Color ; in.ps.Inv.Color; in.ps.Con.Color];
colormap(h, c_map);

xt = unique(M.(PAR_VARX));
yt = unique(M.(PAR_VARY));
dx = xt(2) - xt(1);
dy = yt(2) - yt(1);
[xMesh, yMesh] = meshgrid(xt, yt);
xPatch = [xMesh(:)' - dx/2 ; xMesh(:)' + dx/2 ; xMesh(:)' + dx/2 ; xMesh(:)' - dx/2]; % column is a point, row is a vertix of square centred at point, going anticlockwise from SW corner
yPatch = [yMesh(:)' - dy/2 ; yMesh(:)' - dy/2 ; yMesh(:)' + dy/2 ; yMesh(:)' + dy/2];

M_actions = table2array(M(:, in.action_stems));
M.MaxAction = strings(height(M),1); % type of action that is max
M.MaxActionColorValue = zeros(height(M),1); % color value for max action according to cmap
for i = 1:height(M_actions) % add dots overlay   
    action_max = max(M_actions(i,:));
    stem = in.action_stems{find(M_actions(i,:) == action_max, 1)}; % find dominant action
    M.MaxAction(i) = stem;
    % M.MaxActionColorValue(i) = find(ismember(in.stems, stem)) - 1 + action_max; % color intensity varies with magnitude
    if contains(stem, 'Con')
        M.MaxActionColorValue(i) = 1;
    elseif contains(stem, 'Inv')
        M.MaxActionColorValue(i) = 0;
    elseif contains(stem, 'Rad')
        M.MaxActionColorValue(i) = -1;
    else
        M.MaxActionColorValue(i) = NaN;
        display('Color Error');
    end
end

patch(xPatch, yPatch, M.MaxActionColorValue, 'EdgeColor', 'black'); %, 'LineWidth', 1);
%1 - (1 - in.LIGHTNESS_MAX * (linspace(1,0,in.N_COLOURS))') * (1 - c_max); % range of base color with lightness varying from base color to max lightness parameter
for i = 1:height(M_actions) % add dots overlay    
    stem_ID = extractBefore(M.MaxAction(i),'_');
    ms = in.ps.(stem_ID).MarkerSize;
    if ms > 0
        plot(M.(PAR_VARX)(i), M.(PAR_VARY)(i), '*', 'MarkerSize', ms, 'Color', in.ps.(stem_ID).MarkerEdgeColor);
    end
end

% set(gca, 'CLim', [0 numel(in.stems)]); 
set(gca, 'CLim', [-1 1]);
axis square;
xlim([min(xPatch(:)) max(xPatch(:))]);
ylim([min(yPatch(:)) max(yPatch(:))]);
set(gca, 'FontSize',12);

if strcmpi(in.PLOT_GROUP_BY, 'Individual')
    labelwrap(h, 'Dominant Action (end of exp)', 'sgtitle');
    labelwrap(h, in.plot_subtitle, 'title', 'tex');    
else % grid plot
    title(in.plot_subtitle);
end
if in.INCLUDE_LABELS
    if in.compact_axis_label
        xlabel(PAR_VARX);
        ylabel(PAR_VARY);
    else
        labelwrap(h, in.ps.(PAR_VARX).DisplayName, 'xlabel', 'tex');
        labelwrap(h, in.ps.(PAR_VARY).DisplayName, 'ylabel', 'tex');
    end
    xticks(xt);
    yticks(yt);    
    cb = colorbar('Ticks', -1 + ((1:3)-1)*2/3 + 2/6, 'TickLabels', ["Radical" "Innovative" "Conventional"], 'Location', 'northoutside');
    cb.Label.String = in.ps.(stem_ID).DisplayName;
else                
    xticks([])
    yticks([])
    cb = colorbar('Ticks', []);
end
hold off;
end