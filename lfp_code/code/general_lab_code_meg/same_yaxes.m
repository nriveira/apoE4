function yAx = same_yaxes(varagin)
% function [yAx, xAx] = same_yaxes(figHandle)
%
% PURPOSE:
%   To make the y axis limits identical in all subplots in a figure.
%
% INPUT:
%   figHandle (optional) = handle for figure to edit. If not included, code
%       will edit most recently created figure.
%
% OUTPUT:
%   yAx = y-axis limits for all subplots.
%
% MMD
% 12/2021
% Colgin Lab

if nargin == 0
    figHandle = gcf;
end %use input or not

% numSub = numel(figHandle.Children); %get number of subplots
findAx = findall(gcf, 'type', 'axes');
numSub = numel(findAx);
for n = 1:numSub
    pos1(n) = round(findAx(n).Position(1),2); %rows - we round bc otherwise, I was someitmes getting a miute differece (e-17) w/ dif length titles
    pos2(n) = round(findAx(n).Position(2),2); %columns
end %all subplots

nCols = numel(unique(pos1));
nRows = numel(unique(pos2));

yMin = inf; %initialize
yMax = -inf;

for sp = 1:numSub
    subplot(nRows, nCols, sp)
    ax = gca;
    
    if yMax < ax.YLim(2)
        yMax = ax.YLim(2);
    end %y max
    if yMin > ax.YLim(1)
        yMin = ax.YLim(1);
    end %y min
    
end %subplots to check  min and max

for sp = 1:numSub
    subplot(nRows, nCols, sp)
    ylim([yMin yMax])
end %subplot to change axes

yAx = [yMin yMax];

end %function