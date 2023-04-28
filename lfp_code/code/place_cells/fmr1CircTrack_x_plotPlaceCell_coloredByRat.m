function fmr1CircTrack_x_plotPlaceCell_coloredByRat(group)
% function fmr1CircTrack_x_plotPlaceCell_coloredByRat(group)
% 
% PURPOSE:
%   Plots the place cell properties with different colors for each rat in
%   dotplots.
% 
% INPUT:
%   group = data struct
% 
% OUTPUT:
%   Figures.
% 
% MMD
% Colgin Lab
% 04/2022

%% OPTIONS

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\placeCellProperties\coloredByRat';
saveOrNot = 1;

cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};

spatBinSz = 4;

minFr = 1; %for units to be included

%% INITIALIZE

spatInfo = cell(3,2); %by group
sparsity = cell(3,2);
meanFirRate = cell(3,2);
peakFieldFirRate = cell(3,2);
meanFieldFirRate = cell(3,2);
fieldSize = cell(3,2);
fieldPerCell = cell(3,2);

cols = cell(3,2);
cols{1,1} = 'Blue';
cols{2,1} = 'LightBlue';
cols{3,1} = 'MidnightBlue';
cols{1,2} = 'Red';
cols{2,2} = 'Salmon';
cols{3,2} = 'DarkRed';

convFact = (2*pi*50)/(360/spatBinSz); %cm per bin, because track is 100 cm in diameter and ratemap has 72 bins

jitter = 0.1; %for dotplots

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap) < minFr
                    continue %to next unit
                end %doens't reach max firing rate
                    rateMap = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap; %get AVG ratemap across all begins
                    timePerBin = group(g).rat(r).day(d).xAllBeginTpb; %same, all begins
                    
                    tmpSpatInfo = get_spatial_info(rateMap, timePerBin);
                    spatInfo{r,g} = [spatInfo{r,g} tmpSpatInfo];
                    
                    tmpSparse = get_spatial_sparsity(rateMap, timePerBin);
                    sparsity{r,g} = [sparsity{r,g} tmpSparse];
                    
                    if isempty(group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf) %only cells with place fields from avg ratemap
                        continue %to next cell
                    end %no pf
                    
                    fieldPerCell{r,g} = [fieldPerCell{r,g} length(group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf)];
                    
                    tmpPeak = [];
                    tmpFieldMean = [];
                    tmpSize = [];
                    
                    for pf = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf) %if there are multiple place fields, get mean of measures for this unit
                        tmpPeak = [tmpPeak group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf(pf).pkFr];
                        
                        pfInds = group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf(pf).inds;
                        pullBins = rateMap(pfInds);
                        tmpFieldMean = [tmpFieldMean mean(pullBins)];
                        tmpSize = [tmpSize length(pfInds)*convFact];
                    end %pf
                    
                    peakFieldFirRate{r,g} = [peakFieldFirRate{r,g} mean(tmpPeak)];
                    meanFieldFirRate{r,g} = [meanFieldFirRate{r,g} mean(tmpFieldMean)];
                    
                    meanFirRate{r,g} = [meanFirRate{r,g} mean(rateMap)];
                    fieldSize{r,g} = [fieldSize{r,g} mean(tmpSize)];
            end %unit
        end %day
    end %rat
end %group

cd(saveDir)
keyboard

%% SPAT INFO

figtitle = 'SpatialInformation';
figure('Name', figtitle)

hold on;
h = []; %line handle
leg = {}; %legend
for g = 1:2 
    for r = [3 2 1]
        if ~isempty(spatInfo{r,g})
            tmpH = dotplot(g, spatInfo(r,g), jitter, [rgb(cols{r,g})], [rgb('Black')]);
            h = [h; tmpH];
            tmpLegText = [groupNames{g} ' rat ' num2str(r) ' (n = ' num2str(length(spatInfo{r,g})) ' cells)'];
            leg = [leg; tmpLegText];
        end %there's data
    end %rat
end %group

xlim([0 3])
xticks(1:2)
xticklabels(groupNames)

ylabel('Spatial information (bits/spike)')

for g = 1:2
    tmpMean = mean(horzcat(spatInfo{g,:}));
    line([g-0.25 g+0.25], [tmpMean tmpMean], 'Color', 'Black', 'LineWidth', 3) %mean lines should be longer than the jitter so you can see them
end

legend(h, leg, 'Location', 'Northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% SPAT SPARSITY

figtitle = 'SpatialSparsity';
figure('Name', figtitle)

hold on;
h = []; %line handle
leg = {}; %legend
for g = 1:2 
    for r = [3 2 1]
        if ~isempty(sparsity{r,g})
            tmpH = dotplot(g, sparsity(r,g), jitter, [rgb(cols{r,g})], [rgb('Black')]);
            h = [h; tmpH];
            tmpLegText = [groupNames{g} ' rat ' num2str(r) ' (n = ' num2str(length(sparsity{r,g})) ' cells)'];
            leg = [leg; tmpLegText];
        end %there's data
    end %rat
end %group

xlim([0 3])
xticks(1:2)
xticklabels(groupNames)

ylim([0 1])
ylabel('Spatial sparsity')

for g = 1:2
    tmpMean = mean(horzcat(sparsity{g,:}));
    line([g-0.25 g+0.25], [tmpMean tmpMean], 'Color', 'Black', 'LineWidth', 3) %mean lines should be longer than the jitter so you can see them
end

legend(h, leg, 'Location', 'Northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIELD SIZE

figtitle = 'PlaceFieldSize';
figure('Name', figtitle)

hold on;
h = []; %line handle
leg = {}; %legend
for g = 1:2 
    for r = [3 2 1]
        if ~isempty(fieldSize{r,g})
            tmpH = dotplot(g, fieldSize(r,g), jitter, [rgb(cols{r,g})], [rgb('Black')]);
            h = [h; tmpH];
            tmpLegText = [groupNames{g} ' rat ' num2str(r) ' (n = ' num2str(length(fieldSize{r,g})) ' cells)'];
            leg = [leg; tmpLegText];
        end %there's data
    end %rat
end %group

xlim([0 3])
xticks(1:2)
xticklabels(groupNames)

ylabel('Place field size (cm)')
ylim([0 100])

for g = 1:2
    tmpMean = mean(horzcat(fieldSize{g,:}));
    line([g-0.25 g+0.25], [tmpMean tmpMean], 'Color', 'Black', 'LineWidth', 3) %mean lines should be longer than the jitter so you can see them
end

legend(h, leg, 'Location', 'Northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIELDS PER CELL

figtitle = 'PlaceFieldPerCell';
figure('Name', figtitle)

hold on;
h = []; %line handle
leg = {}; %legend
for g = 1:2 
    for r = [3 2 1]
        if ~isempty(fieldPerCell{r,g})
            tmpH = dotplot(g, fieldPerCell(r,g), jitter, [rgb(cols{r,g})], [rgb('Black')]);
            h = [h; tmpH];
            tmpLegText = [groupNames{g} ' rat ' num2str(r) ' (n = ' num2str(length(fieldPerCell{r,g})) ' cells)'];
            leg = [leg; tmpLegText];
        end %there's data
    end %rat
end %group

xlim([0 3])
xticks(1:2)
xticklabels(groupNames)

ylabel('Place fields/cell')

for g = 1:2
    tmpMean = mean(horzcat(fieldPerCell{g,:}));
    line([g-0.25 g+0.25], [tmpMean tmpMean], 'Color', 'Black', 'LineWidth', 3) %mean lines should be longer than the jitter so you can see them
end

legend(h, leg, 'Location', 'Northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% PEAK FIRING RATE

figtitle = 'PeakPlaceFieldFiringRate';
figure('Name', figtitle)

hold on;
h = []; %line handle
leg = {}; %legend
for g = 1:2 
    for r = [3 2 1]
        if ~isempty(peakFieldFirRate{r,g})
            tmpH = dotplot(g, peakFieldFirRate(r,g), jitter, [rgb(cols{r,g})], [rgb('Black')]);
            h = [h; tmpH];
            tmpLegText = [groupNames{g} ' rat ' num2str(r) ' (n = ' num2str(length(peakFieldFirRate{r,g})) ' cells)'];
            leg = [leg; tmpLegText];
        end %there's data
    end %rat
end %group

xlim([0 3])
xticks(1:2)
xticklabels(groupNames)

ylabel('Peak firing rate (Hz)')

for g = 1:2
    tmpMean = mean(horzcat(peakFieldFirRate{g,:}));
    line([g-0.25 g+0.25], [tmpMean tmpMean], 'Color', 'Black', 'LineWidth', 3) %mean lines should be longer than the jitter so you can see them
end

legend(h, leg, 'Location', 'Northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% PEAK FIRING RATE

figtitle = 'AvgPlaceFieldFiringRate';
figure('Name', figtitle)

hold on;
h = []; %line handle
leg = {}; %legend
for g = 1:2 
    for r = [3 2 1]
        if ~isempty(meanFieldFirRate{r,g})
            tmpH = dotplot(g, meanFieldFirRate(r,g), jitter, [rgb(cols{r,g})], [rgb('Black')]);
            h = [h; tmpH];
            tmpLegText = [groupNames{g} ' rat ' num2str(r) ' (n = ' num2str(length(meanFieldFirRate{r,g})) ' cells)'];
            leg = [leg; tmpLegText];
        end %there's data
    end %rat
end %group

xlim([0 3])
xticks(1:2)
xticklabels(groupNames)

ylabel('Average in-field firing rate (Hz)')

for g = 1:2
    tmpMean = mean(horzcat(meanFieldFirRate{g,:}));
    line([g-0.25 g+0.25], [tmpMean tmpMean], 'Color', 'Black', 'LineWidth', 3) %mean lines should be longer than the jitter so you can see them
end

legend(h, leg, 'Location', 'Northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end


end %function