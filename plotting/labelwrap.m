function labelwrap(h, label, location, interpreter) 
%%% Automatically wrap titles and labels to fit figure width
    if ~exist('interpreter', 'var') || isempty(interpreter)
        interpreter = 'none';
    end

    hUnits_Old = h.Units;
    h.Units = 'characters';
    char_limit = [floor(h.Position(3)/2) floor(h.Position(4)*1.4)]; % character units does not give accurate character dimensions
    label_s = string(label);

    switch location
        case 'sgtitle'            
            sgtitle(LaTeXFormat(textwrap(label_s, char_limit(1))), 'Interpreter', interpreter);
        case 'title'            
            title(LaTeXFormat(textwrap(label_s, char_limit(1))), 'Interpreter', interpreter);
        case 'subtitle'            
            subtitle(LaTeXFormat(textwrap(label_s, char_limit(1))), 'Interpreter', interpreter);
        case 'xlabel'
            xlabel(LaTeXFormat(textwrap(label_s, char_limit(1))), 'Interpreter', interpreter);
        case 'ylabel'
            ylabel(LaTeXFormat(textwrap(label_s, char_limit(2))), 'Interpreter', interpreter);
        case 'zlabel'
            zlabel(LaTeXFormat(textwrap(label_s, char_limit(2))), 'Interpreter', interpreter);
    end
    h.Units = hUnits_Old;
end

function char_cell = LaTeXFormat(char_cell)
    continuing_latex = false;
    for i = 1:numel(char_cell)
        s = char_cell{i};
        if continuing_latex
            s = ['$' s];
        end        
        if mod(sum(s == '$'), 2) == 1
            s = [s '$'];
            continuing_latex = true;
        else
            continuing_latex = false;
        end
        char_cell{i} = s;
    end
end