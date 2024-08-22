function print_as_is(figure_set, fname, format) 
%%% Print figure as it is displayed; i.e. without white space
    if ~exist('figure_set', 'var') || isempty(figure_set) || strcmpi(figure_set, 'all')
        figure_set = findobj('type','figure');
    elseif isnumeric(figure_set)
        figure_set_full = findobj('type','figure');
        figure_set = figure_set_full(ismember([figure_set_full.Number], figure_set)); % figures may be order not according to number
    end

    if ~exist('format', 'var') || isempty(format)
        format = '-dpdf';
        fext = '.pdf';
    else
        fext = '';
    end

    for fi = 1:numel(figure_set)
        h = figure_set(fi);
        if ~exist('fname', 'var') || isempty(fname)
            if isempty(h.Name)
                fname = append("Figure", string(h.Number), fext);
            else
                fname = append(h.Name, fext);
            end
        end        
    
        hUnits_Old = h.Units;
        set(h, 'Units', 'centimeters', 'PaperUnits', 'centimeters');
        hPos = h.Position;
        set(h, 'PaperSize', hPos(3:4), 'Units', hUnits_Old);
        switch format
            case '-dpdf'
                print(h, fname, format, '-fillpage');
            case '-dpng'
                print(h, fname, format,'-r300');
            otherwise
                print(h, fname, format);
        end
        fname = []; % resetting for next figure
    end
end