function fmr1CircTrack_x_make2DRateMaps(group)
% function fmr1CircTrack_x_makeRateMaps(group)
% 
% PURPOSE:
%   Plot the 2D ratemaps for all cells in the FMR1 circle track project.
%   For ilustrative purposes only; all functions use 1D ratemaps.
% 
% INPUT:
%   group = data struct
% 
% OUTPUT:
%   Fn-n: Figures with normalized ratemaps across the 4 begin sessions for
%   each unit. 
% 
% OPTIONS:
%   maxCell = max number of cells for each figure.
% 
% MMD
% 2/2022
% Colgin Lab

%% OPTIONS

saveOrNot = 1; % Set to 1 to save plots in saveDir (1 fig for each group/rat/day)

maxCell = 6; %max number of cells per fig (for legibility)

spatBinSz = 4; %cm

%% INITIALIZE

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\2DrateMaps';

xBnds = [-60 60];
yBnds = [-60 60];

plotOrNot = 0;
velFilt = 1;
durCrit = 1;

numXBins = floor(diff(xBnds)/spatBinSz); % Number of spatial bins is rounded down range of x or y values...
numYBins = floor(diff(yBnds)/spatBinSz); % divided by the spatial bin size

spOpts = 1:4*maxCell; %subplot options
spOpts = reshape(spOpts, 4, maxCell);
spOpts = spOpts';

minVal = -0.01;
maxVal = 1;
cMap = [1 1 1; jet(length(minVal:0.01:maxVal))];

curDir = pwd;
cd(saveDir)

%% GET AND PLOT DATA

for g = 1:2
    fprintf('Group %d\n', g);
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            
            numU = length(group(g).rat(r).day(d).xAllBeginUnitInfo);
            numFigs = ceil(numU/maxCell);
            
            u = 0;
            for f = 1:numFigs
                
                if numFigs == 1
                    figtitle = [group(g).name '_' group(g).rat(r).name '_day' num2str(d)];
                else
                    figtitle = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_' num2str(f) 'of' (num2str(numFigs))];
                end %how many figs
                
                figU = numU-u; %units for this figure
                if figU > maxCell
                    figU = maxCell - figU;
                end
                
                
                figure('Name', figtitle, 'Position', [549 89 923 884])
                hold on;
                
                for figUCntr = 1:maxCell
                    u = u + 1;
                    uMax = 0; %max firing rate
                    if u > numU
                        break
                    end %all units done
                    
                    uID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                    uName = ['TT' num2str(uID(1)) '\_' num2str(uID(2))];
                    allRateMaps = zeros(numXBins, numYBins, 4); %for this unit
                    
                    for b = 1:4
                        spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        
                        [rateMap, spkCnt, timePerBin] = get_2d_ratemap(spkTms, coords, xBnds, yBnds, spatBinSz, plotOrNot, velFilt, durCrit);
                        smRateMap = smooth_2d_ratemap(rateMap);
                        
                        allRateMaps(:,:,b) = smRateMap;
                    end %begin
                    
                    maxFr = max(allRateMaps(:));
                    
                    for b = 1:4
                        subplot(maxCell,4,spOpts(figUCntr,b))
                        tmpMap = allRateMaps(:,:,b) / maxFr;
                        tmpMap(isnan(tmpMap)) = minVal;
                        
                        imagesc(tmpMap)
                        axis xy
                        axis square
                        colormap(cMap)
                        
                        xticks('')
                        yticks('')
                        
                        if b == 1
                            ylabel({[uName]; ['Max = ' num2str(round(maxFr)) ' Hz']; 'Position'});
                        end %if first begin window
                        
                        if figUCntr == maxCell || u == numU %if last unit in this fig
                            xlabel('Position');
                        end
                        
                        if figUCntr == 1
                            title(['Begin ' num2str(b)]);
                        end %if first unit
                        
                    end %begin
                end %units in this fig
                if saveOrNot == 1
                    saveas(gcf, figtitle, 'epsc');
                    saveas(gcf, figtitle, 'fig');
                    saveas(gcf, figtitle, 'png');
                end %save or not
            end %figs
        end %day
    end %rat
end %group



end %function