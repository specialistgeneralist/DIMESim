function res_timeseries_plot(RESFNAME, PAR_CONST_SET, WIDE_FORMAT, MINIMAL_DETAIL)
tic;
    %% Input Extraction
load(RESFNAME, 'inputs', 'res', 'T', 'G', 'NodeActions');
fprintf('Loading finished: %.2f s\n', toc);

if ~exist('PAR_CONST_SET', 'var')
    PAR_CONST_SET = [];     % -- string ("") array of names of constant parameters to be printed in figure. e.g. ["model_variant", "nu", "R"]
end
if ~exist('WIDE_FORMAT', 'var') || isempty(WIDE_FORMAT)
    WIDE_FORMAT = false;
end
if ~exist('MINIMAL_DETAIL', 'var') || isempty(MINIMAL_DETAIL)
    MINIMAL_DETAIL = false;
end

RESFNAME_split = split(RESFNAME, {filesep, '.mat'});
simulation_name = RESFNAME_split{end-1};
n_agents = inputs.n;
rolling_average_window = inputs.rolling_average_window;
network_plot_times = inputs.network_plot_times; 
t_set = 1:T;
n_npt = numel(network_plot_times);

    %% Plot Setup
ps = plot_settings_common('line only');
ps_graph = plot_settings_common();
if WIDE_FORMAT
    figure_size = [720 285];
else
    figure_size = [500 380];
end

t_plot = t_set; %round(exp(linspace(log(1), log(T), 100)));
plot_subtitle = "";
for par = PAR_CONST_SET
    par_value = inputs.(par);
    if iscell(par_value)
        par_value = par_value{1};    % convert from cell form. usually for strings
    end
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
    plot_subtitle = append(plot_subtitle, par_string, '; ');
end
% plot_subtitle = append(plot_subtitle, "$");
if isempty(PAR_CONST_SET)
    timeseries_subtitle = "";
else
    timeseries_subtitle = append(plot_subtitle, 'RA_{window} =', num2str(rolling_average_window));
end
fprintf('Plot setup finished: %.2f s\n', toc);

    %% Plotting
fi = 0;

%{/
fi = fi+1; fh = figure_subplot(fi, [1 1], figure_size); clf
    fh.Name = [simulation_name '_Action'];
    set(gca, 'FontSize', 10.5);
    % set(fh, 'defaultAxesColorOrder', [0 0 0; ps.Ina.Color]); % using inactive color for switching
    hold on;    
    % yyaxis left    
    % for stem_ID = ["Con" "Rad" "InaConLat" "InaRadLat" "C"]
    for stem_ID = ["ConHab" "RadHab" "Inv" "InaConHab" "InaRadHab" "InaInv"]
        stem = stem_ID + "_ra";
        h = plot(t_plot, [res(t_plot).(stem)] * 100);
        set(h, ps.(stem_ID))
    end    
    xlim([1 T]);
    ylim([0, 100]);
    % xscale('log');
    set(gca,'XScale','log');        % -- legacy compatible
    if MINIMAL_DETAIL
        xlabel('Timestep, t');
        ylabel('Population Fraction');
        legend('show', 'Location', 'NorthEast');
        legend('boxoff');
    else % sufficent detail as stand-alone image
        xlabel('Timestep, t')
        labelwrap(fh, 'Rolling average in time of action percentage averaged over agents and replicates', 'ylabel');
        labelwrap(fh, 'Action Rolling Average (Timeseries)', 'title');
        labelwrap(fh, timeseries_subtitle, 'subtitle', 'TeX');
        % labelwrap(fh, sprintf('Showing %d/%d points chosen uniformly on log-scale, with rolling average window of %d' ...
        %     , numel(t_plot), numel(t_set), ROLLING_AVERAGE_WINDOW), 'subtitle');
        legend('show', 'Location', 'NorthEast');
        legend('boxoff'); 
    end
    grid on;       
    hold off;
    set(gca,'FontSize', 14, ...
        'Xtick', 10.^[0:4], ...
        'YTick', 0:10:100)
    set(gcf,'Position',[1   623  392   324])       % -- ensure the same every time
    legend off
    % print as png, with resolution 300
    print([simulation_name '_Action'],'-dpng','-r300');

fi = fi+1; fh = figure_subplot(fi, [1 1], figure_size); clf
    fh.Name = [simulation_name '_DIME']; 
    set(gca, 'FontSize', 10.5);
    hold on;
    for stem_ID = ["D", "I", "M", "E"]        
        stem = stem_ID + "_ra";
        h = plot(t_plot, [res(t_plot).(stem)]);
        set(h, ps.(stem_ID));
    end
    xlim([1 T]);
    ylim([0, 100]);    
    % xscale('log');
    set(gca,'XScale','log');        % -- legacy compatible
    if MINIMAL_DETAIL
        xlabel('Timestep, t');
        ylabel('Mean DIME Values');
        legend('show', 'Location', 'NorthWest');
        legend('boxoff');
    else % sufficent detail as stand-alone image    
        xlabel('Timestep, t');
        labelwrap(fh, 'Rolling average in time of DIME averaged over agents and replicates', 'ylabel');    
        labelwrap(fh, 'DIME Rolling Average (Timeseries)', 'title');
        labelwrap(fh, timeseries_subtitle, 'subtitle', 'TeX');
        % labelwrap(fh, sprintf('Showing %d/%d points chosen uniformly on log-scale, with rolling average window of %d' ...
            % , numel(t_plot), numel(t_set), ROLLING_AVERAGE_WINDOW), 'subtitle');
        legend('show', 'Location', 'NorthWest');
        legend('boxoff');
    end
    grid on;    
    hold off;
    set(gca,'FontSize', 14, ...
        'Xtick', 10.^[0:4], ...
        'YTick', 0:10:100)
    set(gcf,'Position',[1   623   392   324])       % -- ensure the same every time
    legend off
    % print as png, with resolution 300
    print([simulation_name '_DIME'],'-dpng','-r300');

%}  
fi = fi+1; fh = figure_subplot(fi, [1 n_npt], figure_size); clf;
    fh.Name = [simulation_name '_Network'];
    % set(fh, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.4, 0.6]);
    set(fh, 'Renderer','painters','Color',[1 1 1])
    set(axes,'YColor','none','XColor','none')    
    for i = 1: n_npt
        t = network_plot_times(i);
        h = subplot(1, n_npt, i);  
        hold on;
        set(h, 'Color',[1 1 1])        
        set(h,'YColor','none','XColor','none')
        NodeColors = zeros(n_agents, 3);
        NodeMarkers = strings(n_agents, 1);
        for j =1:n_agents
            NodeAction = NodeActions(i,j);
            NodeColors(j,:) = ps_graph.(NodeAction).Color;
            NodeMarkers(j) = string(ps_graph.(NodeAction).Marker);
        end        
        plot(G, 'NodeColor', NodeColors, 'Marker', NodeMarkers, ...
            'LayOut', 'force', 'NodeLabel', [], ...
            'MarkerSize', 10, 'EdgeColor', 'Black', 'EdgeAlpha', 0.3);   %'LineWidth', 1.5
        % if ~MINIMAL_DETAIL
            title(sprintf('Time, $t = %d$', t), 'Interpreter', 'LaTeX');        
        % end
    end
    if ~MINIMAL_DETAIL
        sgtitle({append('Snapshots of Action Network at Different Times: ', plot_subtitle), ...
            'blue = Conventional, teal = Innovator, green = Radical, o = Active, * = Latent'}, ...
            'interpreter', 'tex')
        % sgtitle({append('Snapshots of Action Network at Different Times: ', plot_subtitle), ...
        %     '{\color{blue}o} = Active Conventional, {\color{blue}*} = Latent Conventional, {\color{green}o} = Active Radical, {\color{green}*} = Latent Radical'}, ...
        %     'interpreter', 'tex')
    end
    % labelwrap(fh, plot_subtitle, 'subtitle', 'LaTeX');

fprintf('Plot finished: %.2f s\n', toc);