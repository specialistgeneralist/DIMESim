function f = figure_subplot(flabel, plot_rc, ideal_size)
%%% Adjust the position and size of figure window for subplot grid so that
%%% each subplot is in good size. Clears content if existing.
%%% ideal_size of full figure in pixels

    if ismember(class(flabel), {'matlab.ui.figure', 'matlab.graphics.axis.Axes'}) % option to return subplot handles
        f = flabel;
        return
    end

    figure_set = findobj('type','figure');    
    f = figure(flabel);
    if ~isempty(figure_set) && ismember(flabel, [figure_set.Number]) % if exists, no need to adjust size and position
        clf(f);
        return
    end
    if ~exist('plot_rc', 'var') || isempty(plot_rc)
        plot_rc = [1 1];
    end
    if ~exist('ideal_size', 'var') || isempty(ideal_size)
        ideal_length = 350; % pixels. for one subplot    
        ideal_size = ideal_length * plot_rc(end:-1:1);
    end
    
    screen_size = get(0,'ScreenSize');    
    scaling_ratio = ideal_size / max(ideal_size);
    if ideal_size(1) > 0.9 * screen_size(3)
        if ideal_size(2) > 0.9 * screen_size(4)
            if ideal_size(1) / 0.9 * screen_size(3) > ideal_size(2) / 0.9 * screen_size(4)
                figure_size = 0.9 * screen_size(3) * scaling_ratio;
            else
                figure_size = 0.9 * screen_size(4) * scaling_ratio;
            end
        else
            figure_size = 0.9 * screen_size(3) * scaling_ratio;
        end
    elseif ideal_size(2) > 0.9 * screen_size(4)
        figure_size = 0.9 * screen_size(4) * scaling_ratio;
    else
        figure_size = ideal_size;
    end
    
    screen_centre = [mean(screen_size([1, 3])), mean(screen_size([2, 4]))];
    fUnitsOld = f.Units;
    f.Units = 'Pixels';
    f.Position = [(screen_centre - figure_size/2) figure_size];
    f.Units = fUnitsOld;
end