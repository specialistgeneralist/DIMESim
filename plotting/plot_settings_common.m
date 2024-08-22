function [ps, AdditionalCSet] = plot_settings_common(mode, additional_ncolors)

if ~exist('mode', 'var')
    mode = 'normal';
end
if ~exist('additional_ncolors', 'var')
    additional_ncolors = 0;
end

ps = struct();
ps.default = struct('DisplayName', 'none', 'LineWidth', 2, 'LineStyle', '-'...
    , 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'auto');
switch lower(mode)
    case 'line only'
        ps.default.Marker = 'None';    
    case 'bar'
        ps.default.LineStyle = '-';
        ps.default.LineWidth = 2;
        ps.default.EdgeColor = 'k';
        ps.default.HatchStyle = 'fill';
        ps.default.HatchAngle = 45;
        ps.default.HatchColor = 'auto'; % edge color
    case 'contourf'
        ps.default.MarkerSize = 0;
        ps.default.LineWidth = 1.5;
        ps.default.EdgeColor = 'k';
        ps.default.HatchStyle = 'fill';
        ps.default.HatchAngle = 45;
        ps.default.HatchColor = 'k';
        ps.default.HatchLineStyle = '-';
end

c_lines = lines(2);
colors_excluded = [1 1 1 ; 0 0 0; 1 0 0]; % distinguish from white, black, red
% cset = distinguishable_colors(8, colors_excluded);
cset = [0 70 255; ... % conventional
    0 160 0 ; ... % radical
    0 200 200 ; ... % innovating    
    c_lines(2,:) * 255 ; ... % D
    240 230 0 ; ... % I
    ]/255; % M, E
cset = [cset ; distinguishable_colors(7 - size(cset,1), [colors_excluded ; cset])];
AdditionalCSet = distinguishable_colors(additional_ncolors, [colors_excluded ; cset]);

% MATLAB copies struct by value, not reference. So new copy is made for
% each    
% auto edge color copies line color
ps.Con = assignfield(ps.default, 'DisplayName', 'Conventional', 'Color', cset(1,:), 'MarkerFaceColor', cset(1,:));
ps.ConHab = assignfield(ps.Con, 'DisplayName', 'Active Conventional Protestors'); % Consistent
ps.ConInv = assignfield(ps.Con, 'DisplayName', 'Innovating Conventional');
ps.InaConLat = assignfield(ps.Con, 'DisplayName', 'Latent Conventional', 'LineStyle', '--', 'Marker', '*', 'MarkerEdgeColor', 'white', 'MarkerFaceColor', 'none');
ps.InaConHab = assignfield(ps.InaConLat, 'DisplayName', 'Latent Conventional Protestors'); % Consistent
ps.InaConInv = assignfield(ps.InaConLat, 'DisplayName', 'Innovating Latent Conventional');

ps.Rad = assignfield(ps.default, 'DisplayName', 'Radical', 'Color', cset(2,:), 'MarkerFaceColor', cset(2,:));
ps.RadHab = assignfield(ps.Rad, 'DisplayName', 'Active Radical Protestors'); % Consistent 
ps.RadInv = assignfield(ps.Rad, 'DisplayName', 'Innovating Radical');
ps.InaRadLat = assignfield(ps.Rad, 'DisplayName', 'Latent Radical', 'LineStyle', '--', 'Marker', '*', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'none'); 
ps.InaRadHab = assignfield(ps.InaRadLat, 'DisplayName', 'Latent Radical Protestors'); % Consistent
ps.InaRadInv = assignfield(ps.InaRadLat, 'DisplayName', 'Innovating Latent Radical');

ps.Inv = assignfield(ps.default, 'DisplayName', 'Active Innovators', 'Color', cset(3,:), 'MarkerFaceColor', cset(3,:));
ps.InaInv = assignfield(ps.Inv, 'DisplayName', 'Latent Innovators', 'LineStyle', '--', 'Marker', '*', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'none');
% ps.C = assignfield(ps.default, 'DisplayName', 'Innovating', 'Color', cset(8,:), 'MarkerFaceColor', cset(8,:));
% ps.Ina = assignfield(ps.default, 'DisplayName', 'Inactive', 'Color', cset(3,:), 'MarkerFaceColor', cset(3,:));

ps.DominantAction = assignfield(ps.default, 'DisplayName', 'Dominant Action');
ps.DominantAction3D = assignfield(ps.default, 'DisplayName', 'Dominant Action');

ps.D = assignfield(ps.default, 'DisplayName', 'Disidentification', 'Color', cset(4,:), 'MarkerFaceColor', cset(4,:));
ps.I = assignfield(ps.default, 'DisplayName', 'Innovation', 'Color', cset(5,:), 'MarkerFaceColor', cset(5,:));
ps.M = assignfield(ps.default, 'DisplayName', 'Moralisation', 'Color', cset(6,:), 'MarkerFaceColor', cset(6,:));
ps.E = assignfield(ps.default, 'DisplayName', 'Energisation', 'Color', cset(7,:), 'MarkerFaceColor', cset(7,:));

switch lower(mode)
    case 'line only'        
        for prop = fields(ps)'
            ps.(prop{1}).Marker = 'None';
        end
    case 'bar'        
        ps.ConHab = assignfield(ps.ConHab, 'HatchStyle', 'fill', 'HatchAngle', 0, 'HatchColor', [0 0 0]);
        ps.RadHab = assignfield(ps.RadHab, 'HatchStyle', 'fill', 'HatchAngle', 0, 'HatchColor', [0 0 0]);
        ps.Inv = assignfield(ps.Inv, 'HatchStyle', 'fill', 'HatchAngle', 0, 'HatchColor', [0 0 0]);
        % ps.ConInv = assignfield(ps.ConInv, 'HatchStyle', 'cross', 'HatchColor', ps.ConHab.HatchColor);       
        % ps.RadInv = assignfield(ps.RadInv, 'HatchStyle', 'cross', 'HatchColor', ps.RadHab.HatchColor);
        ps.InaConHab = assignfield(ps.InaConHab, 'HatchStyle', 'single', 'HatchAngle', 0, 'HatchColor', ps.ConHab.HatchColor);
        ps.InaRadHab = assignfield(ps.InaRadHab, 'HatchStyle', 'single', 'HatchAngle', 0, 'HatchColor', ps.RadHab.HatchColor);
        ps.InaInv = assignfield(ps.InaInv,  'HatchStyle', 'single', 'HatchAngle', 0, 'HatchColor', ps.Inv.HatchColor);
        % ps.InaConLat = assignfield(ps.InaConLat, 'EdgeColor', ps.Ina.Color, 'HatchStyle', 'speckle', 'HatchColor', ps.ConHab.HatchColor);        
        % ps.InaConInv = assignfield(ps.InaConInv, 'EdgeColor', ps.Ina.Color, 'HatchStyle', 'cross', 'HatchColor', ps.ConHab.HatchColor);
        % ps.InaRadLat = assignfield(ps.InaRadLat, 'EdgeColor', ps.Ina.Color, 'HatchStyle', 'speckle', 'HatchColor', ps.RadHab.HatchColor);        
        % ps.InaRadInv = assignfield(ps.InaRadInv, 'EdgeColor', ps.Ina.Color, 'HatchStyle', 'cross', 'HatchColor', ps.RadHab.HatchColor);
    case 'contourf'
        ms_dots = 5;
        ps.InaConLat.MarkerSize = ms_dots;
        ps.InaRadLat.MarkerSize = ms_dots;         
        ps.InaConHab = assignfield(ps.InaConHab, 'HatchStyle', 'single', 'HatchLineStyle', ':');
        ps.InaRadHab = assignfield(ps.InaRadHab, 'HatchStyle', 'single', 'HatchLineStyle', ':');
        ps.InaInv = assignfield(ps.InaInv, 'HatchStyle', 'single', 'HatchLineStyle', ':');
end

ps.p = struct('DisplayName', 'Probability of failure signal from authorities, p', 'SubtitleName', 'p');
ps.F = struct('DisplayName', 'Threshold for individual reframing, F', 'SubtitleName', 'F');
ps.nu = struct('DisplayName', 'Threshold for collective reframing, \phi', 'SubtitleName', '\phi');
ps.R = struct('DisplayName', 'Length of free communication, R', 'SubtitleName', 'R');
ps.initial_action = struct('DisplayName', 'Initial Action', 'SubtitleName', 'x(0)');
ps.pFnuR = struct('DisplayName', 'Parameter Scenarios', 'SubtitleName', 'pFnuR');
ps.n_replicates = struct('SubtitleName', 'n_{replicates}');

end


function s = assignfield(s, varargin)
% setfield can only assign a single field with multiple indexing levels

for i = 1:2:numel(varargin)
    s.(varargin{i}) = varargin{i+1};
end

end

