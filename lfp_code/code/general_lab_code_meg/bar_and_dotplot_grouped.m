function bar_and_dotplot_grouped(data, yLab, xLab, group1Names, group2Names, cols)
% function bar_and_dotplot_grouped(data, yLab, xLab, group1Names, group2Names, cols)
%
% PURPOSE:
%   To make bar and dot plots showing the same data across two subplots.
%
% INPUTS:
%   data = cell array with individual data points for each group, with
%       dimensions group1 (different colors, grouped together on x-axis) x
%       group2
%   yLab = text string for the ylabel (data measure and units)
%       ex: yLab = 'Spatial information (bits/spike)';
%   xLab = text string for the x label (of which group2Names correspond)
%       ex: xLab = 'Event type';
%   group1Names = cell array with names of each group that is grouped
%       together on x-axis (dim 1 of data)
%       ex: group1Names = {'WT', 'FXS'};
%   group2Names = cell array with names of each group that is separated
%       across x-axis (dim 2 of data)
%       group2Names = {'Forward', 'Reverse'};
%   cols = cell array with names of colors for each of group 1
%       ex: cols = {'Blue', 'Red'};
%
% OUTPUT:
%   Plotted data, with bar plot as subplot 1 and dotplot as subplot 2.
%
% MMD
% 12/2021
% Colgin Lab

%% CHECK

if min(size(data)) == 1
    error('Use bar_and_dotplot function for non-grouped data')
end

if size(data,1) ~= length(group1Names) || size(data,2) ~= length(group2Names) || size(data,1) ~= length(cols)
   error('Group names and colors must match data size') 
end

%% INITIALIZE

jitter = 0.1;  %for dotplot

colMatrix = []; %initialize
for c = 1:length(cols)
    colMatrix = [colMatrix; rgb(cols{c})]; %fill in for dotplot
end %cols

avgData = cellfun(@mean, data)'; %mean
semData = cellfun(@semfunct, data)'; %sem

%% PLOT

subplot(1,2,1)

bgraph = bar(avgData, 'FaceColor', 'Flat');
hold on;

ngroups = size(data, 2);
nbars = size(data, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));

for g = 1:size(data,1)
    bgraph(g).CData = rgb(cols{g});
    
    x = (1:ngroups) - groupwidth/2 + (2*g-1) * groupwidth / (2*nbars);
    er = errorbar(x, avgData(:,g), semData(:,g));
    er.Color = [0 0 0];
    er.LineStyle = 'None';
end %group

ylabel(yLab)
xticklabels(group2Names)
xlabel(xLab)
legend(group1Names, 'Location', 'northeastoutside')

subplot(1,2,2)
lh = zeros(1,2);
for g = 1:size(data,1)
    
    x = (1:ngroups) - groupwidth/2 + (2*g-1) * groupwidth / (2*nbars);
    h = dotplot(x, data(g,:), jitter, repmat(rgb(cols{g}), length(data(g,:)), 1), repmat(rgb('Black'), length(data(g,:)), 1));
    lh(g) = h(1);
    
    for g2 = 1:size(data,2)
         line([x(g2)-jitter x(g2)+jitter], [avgData(g2,g) avgData(g2,g)], 'Color', 'Black', 'LineWidth', 3)
    end %group 2
end %group

ylabel(yLab)
xticklabels(group2Names)
xlabel(xLab)
legend(lh, group1Names, 'Location', 'northeastoutside')

end %function