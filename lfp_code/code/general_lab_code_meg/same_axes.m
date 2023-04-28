function [yAx, xAx] = same_axes(varargin)
% function [yAx, xAx] = same_axes(sps, figHandle)
% 
% PURPOSE:
%   To make the x and y axis limits identical in all subplots in a figure.
% 
% INPUT:
%   sps (optional) = which subplots to check across; default is all
%   figHandle (optional) = handle for figure to edit. If not included, code
%       will edit most recently created figure.
%
% OUTPUT:
%   yAx = y-axis limits for all subplots.
%   xAx = x-axis limits for all subplots.
% 
% MMD
% 12/2021
% Colgin Lab

if nargin < 2
    figHandle = gcf;
end %use input or not

% numSub = numel(figHandle.Children); %get number of subplots
findAx = findall(gcf, 'type', 'axes');
numSub = numel(findAx);

if nargin == 0
    sps = 1:numSub;
else
    sps = varargin{1};
end %use input or not

for n = 1:numSub
    pos1(n) = round(findAx(n).Position(1),2); %rows - we round bc otherwise, I was someitmes getting a miute differece (e-17) w/ dif length titles
    pos2(n) = round(findAx(n).Position(2),2); %columns 
end %all subplots

nCols = numel(unique(pos1));
nRows = numel(unique(pos2));

yMin = inf; %initialize
yMax = -inf;

xMin = inf; %initialize
xMax = -inf;

for spInd = 1:length(sps)
    sp = sps(spInd);
    subplot(nRows, nCols, sp)
    ax = gca;
    
    if yMax < ax.YLim(2)
        yMax = ax.YLim(2);
    end %y max
    if yMin > ax.YLim(1)
        yMin = ax.YLim(1);
    end %y min
    
    if xMax < ax.XLim(2)
        xMax = ax.XLim(2);
    end %x max
    if xMin > ax.XLim(1)
        xMin = ax.XLim(1);
    end %x min
end %subplots to check  min and max

for spInd = 1:length(sps)
    sp = sps(spInd);
    subplot(nRows, nCols, sp)
    ylim([yMin yMax])
    xlim([xMin xMax])
end %subplot to change axes

yAx = [yMin yMax];
xAx = [xMin xMax];

end %function