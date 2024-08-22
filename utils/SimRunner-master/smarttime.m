function [sfstr] = smarttime(fstr)

% SMARTTIME Provides human-readable finish times, relative to now.
%    [SFSTR] = SMARTTIME(FSTR) converts a final finish time string FSTR into a 
%    'smart' finish time string SFSTR by applying an Apple-like relative time 
%    binning (e.g. 'any time now' 'less than a minute');
%
% See also DATESTR DATENUM DATEVEC NOW

% Author: SA, 16 Sep 2005 (UNSW)

% initial data
now_dn = now;
ft_dn  = datenum(fstr);
diff_dn = ft_dn - now_dn;

% get cut-off points for FT_DN in datenums
diff_10s = datenum([0 0 0 0 0 10]);
diff_60s = datenum([0 0 0 0 1 0 ]);
diff_90s = datenum([0 0 0 0 1 30]);
diff_midnight = datenum('24:00:00') - datenum(datestr(now_dn,13));
diff_tom_mdnght = diff_midnight + datenum([0 0 1 0 0 0]);
diff_week_away = datenum([0 0 7 0 0 0]);

% construct intelligent strings
if diff_dn <= diff_10s
	sfstr = 'any time now';
elseif diff_dn <= diff_60s
	sfstr = 'less than a minute';
elseif diff_dn <= diff_90s
	sfstr = 'about a minute';
elseif diff_dn <= diff_midnight
	sfstr = ['today at ' datestr(ft_dn,16)];
elseif diff_dn <= diff_tom_mdnght
	sfstr = ['tomorrow at ' datestr(ft_dn,16)];
elseif diff_dn <= diff_week_away
	sfstr = ['this coming ' datestr(ft_dn,8) ' at ' datestr(ft_dn,16)];
else
	sfstr = datestr(ft_dn,0);
end

