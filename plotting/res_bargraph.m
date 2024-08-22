function res_bargraph(RESFNAME, PAR_VAR, varargin)

%% Analyse final outcomes after varying `par_varied`.

% -- constant settings
in = struct('RESFNAME', RESFNAME);
ps = plot_settings_common('bar');

% -- input parameters
PAPER_FIGURE = []; % -- (3--5, -1 -- -4) Number of relevant figure from paper. Loads hardcoded settings for all parameters

COMPARTMENTALIZE = false;   % -- whether to move input and output files related to this run into a separate folder
in.DO_CONF_INTS = true;     % -- whether to plot confidence intervals
in.MARK_UNCONVERGED = false;     % -- whether to mark unconverged points
ORIENTATION = 'vertical'; % -- vertical or horizontal bars
PAR_CONST_SET = [];     % -- string ("") array of names of constant parameters to be printed in figure. e.g. ["model_variant", "nu", "R"]
PRINT_PLOT = '';        % -- print format for plots. if '', then don't print

% -- handle varargin as name-value pairs. Name is not case-sensitive
for i = 1:2:length(varargin)
    NAME = lower(varargin{i});
    VALUE = varargin{i+1};
    switch NAME
        case 'compartmentalize'
            COMPARTMENTALIZE = VALUE;
        case 'do_conf_ints'
            in.DO_CONF_INTS = VALUE;
        case 'mark_unconverged'
            in.MARK_UNCONVERGED = VALUE;
        case 'orientation'
            ORIENTATION = VALUE;
        case 'par_const_set'
            PAR_CONST_SET = VALUE; 
        case 'print_plot'
            PRINT_PLOT = VALUE;         
    end
end

% -- load results
RESFNAME_var = who('-file', RESFNAME);
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

RESFNAME_split = split(RESFNAME, {filesep, '.mat'});
simulation_name = RESFNAME_split{end-1};
stopping_epsilon = inputs_collated.stopping_epsilon;
% .. by default, we do confidence intervals if there are multiple replicates
if inputs_collated.n_replicates == 1
    in.do_conf_ints = false;
end

% -- checking that constant parameters are really constant over
% simulations
plot_subtitle = '';
for par = PAR_CONST_SET
    par_value = inputs_collated.(par);
    if numel(par_value) > 1 && ~ischar(par_value)
        print([par ' omitted from title since it varies'])
    else
        if iscell(par_value)
            par_value = par_value{1};    % convert from cell form. usually for strings
        end        
        if isequal(par, 'model_variant')
            plot_subtitle = [plot_subtitle par_value(isstrprop(par_value,'upper')) '; '];     % using abbreviation for model variant
        elseif ismember(par, fields(ps))
            plot_subtitle = sprintf('%s%s=%.2g, ', plot_subtitle, ps.(par).SubtitleName, par_value);
        else
            plot_subtitle = sprintf('%s%s=%.2g, ', plot_subtitle, par, par_value);
        end
    end
end
plot_subtitle = [plot_subtitle '. Hatching: None = active; --- = latent'];

in.ORIENTATION = ORIENTATION;
in.ps = ps;
in.plot_subtitle = plot_subtitle;
in.stopping_epsilon = stopping_epsilon;
filename_prefix = erase(RESFNAME, '.mat');

% -- Plot :: confidence intervals underlay as envelope if we have them, mean as line on top
fh = figure_subplot(1, [], [570 570]); clf
    fh.Name = [simulation_name '_DIME']; 
    in.var_stems = {'D_avg' 'I_avg' 'M_avg' 'E_avg'};
    in.plot_title = 'DIME Value (end of exp)';
    in.ylim = [0 100];
    in.yscaling = 1;
    in.ytickformat = 'auto';
    in.ylabel = 'Mean DIME Values';
    in.legend = {}; % default    
    plot_by_par(in, S, M, PAR_VAR, 'grouped')
    switch string(PRINT_PLOT)
        case "png"
            print(1, '-dpng', '-r100', [filename_prefix '_dime.png'])
        case "fig"
            savefig(1, [filename_prefix '_dime.fig'])   % .fig file more useful for custom editting
    end

fh = figure_subplot(2, [], [570 570]); clf
    fh.Name = [simulation_name '_Action'];     
    in.var_stems = flip({'ConHab_f', 'Inv_f', 'RadHab_f', 'InaConHab_f', 'InaInv_f', 'InaRadHab_f'}); % barchart stacks bottom to top
    in.plot_title = 'Action Percentage (end of exp)';
    in.ylim = [0 100];
    in.yscaling = 100;
    in.ytickformat = 'auto'; % 'percentage' option currently not working to display %
    in.ylabel = 'Population Fraction';
    in.legend = flip({'Conventional Protestors' 'Innovators' 'Radical Protestors' '' '' ''});
    plot_by_par(in, S, M, PAR_VAR, 'stacked')
    switch string(PRINT_PLOT)
        case "png"
            print(2, '-dpng', '-r100', [filename_prefix '_actions.png'])
        case "fig"
            savefig(2, [filename_prefix '_actions.fig'])   % .fig file more useful for custom editting
    end

if COMPARTMENTALIZE
    fold = filename_prefix;
    mkdir(fold);
    addpath(fold);
    copyfile([filename_prefix '.txt'], fold)
    movefile(RESFNAME, fold)
    movefile([filename_prefix '_dime.fig'], fold)
    movefile([filename_prefix '_actions.fig'], fold)
end
end


% --------------------------
function plot_by_par(in, S, M, par_var, groupstyle)
    % if in.DO_CONF_INTS
    %     add_conf_intervals(S, in.var_stems, par_var, in)
    % end    

    nstem = numel(in.var_stems);
    nCases = height(M);
    switch par_var
        case 'pFnuR'
            x_ticklabels = strings(nCases, 1);
            for i = 1:nCases
                x_ticklabels(i) = sprintf('p=%.2g; F=%.2g; \\phi=%.2g; R=%d' ...
                    , M.p(i), M.F(i), M.nu(i), M.R(i));
            end
        case 'initial_action'
            x_ticklabels = strings(nCases, 1);
            for i = 1:nCases
                switch M.(par_var){i}
                    case 'all conventional'
                        x_ticklabels(i) = "All Active Conventional";
                    case 'all inactive conventional'
                        x_ticklabels(i) = "All Latent Conventional";
                    case 'all radical'
                        x_ticklabels(i) = "All Active Radical";
                    case 'all inactive radical'
                        x_ticklabels(i) = "All Latent Radical";
                    case 'random'
                        x_ticklabels(i) = "Random";
                end
            end
        otherwise
            x_ticklabels = string(M.(par_var));
    end
    x_coord = 1:nCases; % for plotting as hatchfill2 requires numeric axes data
    y_mat = zeros(nCases, nstem);
    for i = 1:nstem
        this_var = in.var_stems{i};        
        y_mat(x_coord,i) = M.(this_var) * in.yscaling;
    end                           
    
    switch lower(in.ORIENTATION)
        case 'vertical'
            b = bar(x_coord, y_mat, groupstyle);
            hatchAngleOffset = 0;
        case 'horizontal'
            b = barh(x_coord, y_mat, groupstyle);
            hatchAngleOffset = 90;
    end
    
    hold on;
    for i = 1:nstem
        this_var = in.var_stems{i};
        stem_ID = extractBefore(this_var,'_');
        ps_stem = in.ps.(stem_ID);
        b(i).FaceColor = ps_stem.Color;
        b(i).LineWidth = ps_stem.LineWidth;
        b(i).EdgeColor = ps_stem.EdgeColor;
        b(i).DisplayName = ps_stem.DisplayName;
        hatchfill2(b(i), ps_stem.HatchStyle, 'HatchAngle', ps_stem.HatchAngle + hatchAngleOffset...
            , 'HatchColor', ps_stem.HatchColor, 'HatchLineStyle', ps_stem.LineStyle); % adds hatch entry to legend
    end
    hold off;
    
    lgd = legend('show');
    lgd_String_Corrected = repmat({''}, size(lgd.String)); % start with empty array of character cells    
    if isempty(in.legend)
        lgd_String_Corrected(1:2:end) = lgd.String(1:2:end); % keep only bar entries. cannot directly remove hatch entries    
    else % custom legend
        lgd_String_Corrected(1:2:end) = in.legend;
    end
    legend(lgd_String_Corrected);
    if numel(legend().String) <= 4
        legend('Location', 'northoutside', 'Orientation', 'horizontal');
    else
        legend('Location', 'eastoutside', 'Orientation', 'vertical');
    end
    legend('boxoff');
    
    switch lower(in.ORIENTATION)
        case 'vertical'
            ylim(in.ylim);
            xticks(x_coord);
            xticklabels(x_ticklabels);            
            ytickformat(in.ytickformat)    
            xlabel(in.ps.(par_var).DisplayName)
            ylabel(in.ylabel)
        case 'horizontal' % x and y commands refer to on-screen x and y, not input x and y
            xlim(in.ylim);
            yticks(x_coord);
            yticklabels(x_ticklabels);
            xtickformat(in.ytickformat)    
            ylabel(in.ps.(par_var).DisplayName)
            xlabel(in.ylabel)
    end
    title(in.plot_title)
    labelwrap(gcf, in.plot_subtitle, 'subtitle', 'tex');
    % 'FontSize', 14, 'XTick', 0:0.1:1, ...
        grid on         
end



% ------------------------------------------
function add_conf_intervals(S, var_stems, par_var, in, varargin)
    % -- plot confidence intervals as envelope (patch)
    % C = lines(10);
    xx = [S.(par_var); flipud(S.(par_var))];
    xx_t = transform_xticks(S.(par_var), in.XTICKS_TRANSFORMATION);
    xx_t = [xx_t; flipud(xx_t)];
    for i = 1:numel(var_stems)
        this_var = var_stems{i};
        stem_ID = extractBefore(this_var,'_');
        this_var_ci_low = [this_var '_ci_low'];
        this_var_ci_high = [this_var '_ci_high'];
        yy = [S.(this_var_ci_low); flipud(S.(this_var_ci_high))] * in.yscaling;
        if ~isempty(varargin)
            CL_N = varargin{1};
        else
            CL_N = i;
        end
        patch(xx_t, yy, in.ps.(stem_ID).Color, 'FaceAlpha', 0.2, 'EdgeColor', 'none'...
            , 'HandleVisibility','off'), hold on
        
    end
end
