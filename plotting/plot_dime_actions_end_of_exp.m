function plot_dime_actions_end_of_exp(MATFILE)

load(MATFILE)        % -- gives S, M

RESFNAME = MATFILE;

% -- Plot :: confidence intervals underlay as envelope if we have them, mean as line on top
in.plot_by_F = false;  % -- plot all together
in.RESFNAME = RESFNAME;
in.do_conf_ints = true;  % -- plot confidence intervals as envelope

% -- set up common plotting settings
ps = plot_settings_common('line only');

figure(1),clf

    in.var_stems = {'D_avg' 'I_avg' 'M_avg' 'E_avg'};
    in.title = 'DIME (end of exp)';
    in.ylim = [0 100];
    in.plot_rc = [2 2];
    plot_by_p(in, S, M, ps)
    print(1, '-dpng', '-r150', strrep(RESFNAME, '.mat', '_dime.png'))

figure(2),clf
    in.var_stems = {'ConHab_f' 'RadHab_f' 'Inv_f' 'InaConHab_f' 'InaRadHab_f' 'InaInv_f'};
    in.title = 'Actions (end of exp)';
    in.ylim = [0 1];
    in.plot_rc = [1 3];
    in.plot_by_F = false;  % -- plot all together
    plot_by_p(in, S, M, ps)
    print(2, '-dpng', '-r150', strrep(RESFNAME, '.mat', '_actions.png'))

% ------------------


% --------------------------
function plot_by_p(in, S, M, ps)

% in.plot_by_F = false;
% in.var_stems = {'D_avg' 'I_avg' 'M_avg' 'E_avg'};
% in.do_conf_ints = true;
% in.title = 'DIME (end of exp)';

if in.plot_by_F    % -- plot each DIME var separately
    uF = unique(S.F);
    for k = 1:numel(in.var_stems)
        subplot(in.plot_rc(1),in.plot_rc(2),k)
        this_var = in.var_stems(k);
        for i = 1:numel(uF)
            this_F = uF(i);
            if in.do_conf_ints
                add_conf_intervals(S(S.F==this_F,:), this_var, k)
            end
            add_mean_lines(M(M.F==this_F,:), this_var, k)
        end
        xlabel('p')
        title(sprintf('%s (end of exp)', this_var{1}), 'Interpreter', 'none')
        set(gca,'FontSize', 12,...
            'XTick', 0.05:0.05:0.95, ...
            'YLim', in.ylim)
            axis square
            grid on
        % .. if last subplot, add fname
        if k==numel(in.var_stems)
            add_fname(in.RESFNAME)
        end
    end    
else
    if in.do_conf_ints
        add_conf_intervals(S, in.var_stems)
    end
    add_mean_lines(M, in.var_stems, ps)
    add_fname(in.RESFNAME)
    % legend(in.var_stems, 'Location', 'NorthWest')
    xlabel('p')
    ylabel('Value')
    title(in.title)
    set(gca,'FontSize', 10,...
        'XTick', 0.05:0.05:0.95, ...
        'YLim', in.ylim, ...
        'YTick', in.ylim:(in.ylim(2)-in.ylim(1))/10:in.ylim(2), ...
        'XTickLabelRotation', 45)
        title('')
        axis square
        grid on
end


% -------------------------
function add_fname(RESFNAME)

[pathstr,name,ext] = fileparts(RESFNAME);
text(1.2, 0, [name ext], 'Units', 'normalized',...
    'FontSize', 10,...
    'Rotation', 90,...
    'HorizontalAlignment', 'left',...
    'Interpreter', 'none')


% ------------------------------------------
function add_conf_intervals(S, var_stems, varargin)

% -- plot confidence intervals as envelope (patch)
% C = lines(10);
C = ones(1,3) * 0.8;  % -- light grey
xx = [S.p; flipud(S.p)];
for i = 1:numel(var_stems)
    this_var = var_stems{i};
    this_var_ci_low = [this_var '_ci_low'];
    this_var_ci_high = [this_var '_ci_high'];
    yy = [S.(this_var_ci_low); flipud(S.(this_var_ci_high))];
    if nargin > 2
        CL_N = varargin{1};
    else
        CL_N = i;
    end
    h = patch(xx, yy, C, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
end

% ------------------------------------------
function add_mean_lines(M, var_stems, ps, varargin)

% -- plot means as lines
C = lines(10);
x = M.p;
for i = 1:numel(var_stems)
    this_var = var_stems{i};
    y = M.(this_var);
    if nargin > 3
        CL_N = varargin{1};
        text(x(end)*1.01,y(end), sprintf('F=%0.1f', M.F(end)), 'Color', 'k', 'FontSize', 8)
    else
        CL_N = i;
    end
    h = plot(x, y, '-', 'Color', C(CL_N,:), 'LineWidth', 2);
    this_var = strrep(this_var, '_avg', '');
    this_var = strrep(this_var, '_f', '');
    set(h, ps.(this_var));
end