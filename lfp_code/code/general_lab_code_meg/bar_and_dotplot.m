function bar_and_dotplot(data, yLab, xLabs, cols)
% function bar_and_dotplot(data, yLab, xLabs, cols)
%
% PURPOSE:
%   To make bar and dot plots showing the same data across two subplots.
%
% INPUTS:
%   data = cell array with individual data points for each group
%   yLab = cell with text string for the ylabel (data measure and units)
%       ex: yLab = {'Spatial information (bits/spike)'};
%   xLabs = cell array with xticklabels (names of each group)
%       ex: xLabs = {'WT', 'FXS'};
%   cols = cell array with names of colors for each group
%       ex: cols = {'Blue', 'Red'};
%
% OUTPUT:
%   Plotted data, with bar plot as subplot 1 and dotplot as subplot 2.
%
% MMD
% 12/2021
% Colgin Lab

%% CHECK

if min(size(data)) ~= 1
    error('Use bar_and_dotplot_grouped function for grouped data')
end

if length(data) ~= length(xLabs) || length(data) ~= length(cols) || length(xLabs) ~= length(cols)
   error('Data, x label, and color inputs must be the same size') 
end

%% INITIALIZE

jitter = 0.1;  %for dotplot

colMatrix = []; %initialize
for c = 1:length(cols)
    colMatrix = [colMatrix; rgb(cols{c})]; %fill in for dotplot
end %cols

avgData = cellfun(@mean, data); %mean
semData = cellfun(@semfunct, data); %sem
if isnan(sum(avgData))
    warning('NaNs in dataset')
end

%% PLOT

subplot(1,2,1)

bgraph = bar(avgData, 'FaceColor', 'Flat');
hold on;
errorbar(1:length(data), avgData, semData, 'Color', 'Black', 'LineStyle', 'None')
for g = 1:length(data)
    bgraph.CData(g,:) = rgb(cols{g});
end %group

ylabel(yLab)

set(gca(),'XTickLabel', {sprintf(xLabs{1}) sprintf(xLabs{2})});
%set(gca(),'XTickLabel', {sprintf(xLabs{1})});

subplot(1,2,2)
dotplot(1:length(data), data, jitter, colMatrix, repmat(rgb('Black'), length(data), 1));

for g = 1:length(data)
    line([g-.3 g+.3], [avgData(g) avgData(g)], 'Color', 'Black', 'LineWidth', 3)
end %group

ylabel(yLab)
set(gca(),'XTickLabel', {sprintf(xLabs{1}) sprintf(xLabs{2})});
%set(gca(),'XTickLabel', {sprintf(xLabs{1})});

end %function