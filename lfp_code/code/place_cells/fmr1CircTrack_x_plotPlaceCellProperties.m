function fmr1CircTrack_x_plotPlaceCellProperties(group)
% function fmr1CircTrack_x_plotPlaceCellProp(group)
%
% PURPOSE:
%   Plot place cell firing properties as in Figure 3 of Mably et al. 2016
%   for data collected for the WT/FXS circle track data.
%
% INPUT:
%   group data struct
%
% OUTPUT:
%   F1: Spatial correlation.
%   F2: Rate overlap.
%   F3: Spatial information.
%   F4: Spatial sparsity.
%   F5: Place field size.
%   F6: Place fields per cell.
%   F7: Mean firing rate across begin sessions.
%   F8: Peak in-field firing rate.
%   F9: Average in-field firing rate.
%
% MM Donahue
% 5/2021
% Colgin Lab

%% OPTIONS

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\placeCellProperties';
saveOrNot = 0;

prepForStats = 0;

cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};

spatBinSz = 4;

minFr = 1; %for units to be included

%% INITIALIZE

% Include all CA1 units, as in Mably et al. 2016 (Fig 3)
spatCorr = cell(2,4); %group x combo
rateOverlap = cell(2,4); %same

combos = [1 2; 2 3; 3 4; 1 4]; %combos as listed above
combosTxt = {'1-2', '2-3', '3-4', '1-4'};

% Cells with identified place fields
spatInfo = cell(2,1); %by group
sparsity = cell(2,1);
meanFirRate = cell(2,1);
peakFieldFirRate = cell(2,1);
meanFieldFirRate = cell(2,1);
fieldSize = cell(2,1);
fieldPerCell = cell(2,1);

convFact = (2*pi*50)/(360/spatBinSz); %cm per bin, because track is 100 cm in diameter and ratemap has 72 bins

uCntr = 0; %initialize for stats (ANOVAs)

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap) < minFr
                    continue %to next unit
                end %doens't reach max firing rate
                uCntr = uCntr + 1;
                uRateMaps = zeros(4,360/spatBinSz); %for storing rate maps from all begins for this unit
                
                for b = 1:4
                    uRateMaps(b,:) = group(g).rat(r).day(d).begin(b).unit(u).smRateMap; %get all of the ratemaps for the spat corr comparisons
                end %begin
                
                for c = 1:4
                    b1 = combos(c,1); %begins to compare
                    b2 = combos(c,2);
                    
                    rm1 = uRateMaps(b1,:); %ratemaps to compare
                    rm2 = uRateMaps(b2,:);
                    
                    scMat = corrcoef(rm1,rm2);
                    if ~isnan(scMat(2))
                        spatCorr{g,c} = [spatCorr{g,c} scMat(2)];
                    end
                    
                    fr1 = mean(uRateMaps(b1,:));
                    fr2 = mean(uRateMaps(b2,:));
                    
                    rateRatio = min([fr1 fr2]) / max([fr1 fr2]); % Calculate the rate ratio with the lower FR as numerator
                    if ~isnan(rateRatio)
                        rateOverlap{g,c} = [rateOverlap{g,c} rateRatio];
                    end
                    
                end %combos
                
                if ~isempty(group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf) %only cells with place fields from avg ratemap
                    
                    rateMap = group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap; %get AVG ratemap across all begins
                    timePerBin = group(g).rat(r).day(d).xAllBeginTpb; %same, all begins
                    
                    tmpSpatInfo = get_spatial_info(rateMap, timePerBin);
                    spatInfo{g,1} = [spatInfo{g,1} tmpSpatInfo];
                    
                    tmpSparse = get_spatial_sparsity(rateMap, timePerBin);
                    sparsity{g} = [sparsity{g} tmpSparse];
                    
                    fieldPerCell{g,1} = [fieldPerCell{g,1} length(group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf)];
                    
                    tmpPeak = [];
                    tmpFieldMean = [];
                    tmpSize = [];
                    
                    for pf = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf) %if there are multiple place fields, get mean of measures for this unit
                        tmpPeak = [tmpPeak group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf(pf).pkFr];
                        
                        pfInds = group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf(pf).inds;
                        pullBins = zeros(1,length(pfInds)); %initialize so I know it works right
                        
                        if pfInds(1) < pfInds(end) %it doesn't wrap around 0 degrees
                            pullBins = rateMap(1,pfInds(1):pfInds(end));
                        else %pf does wrap around 0 degrees
                            findWrap = find(pfInds == 1); %find where it crosses 0 degree boundary
                            pullBins(1,1:findWrap-1) = rateMap(1,pfInds(1):pfInds(findWrap-1));
                            pullBins(1,findWrap:end) = rateMap(1,pfInds(findWrap):pfInds(end));
                        end %whether or not place field wraps around 0
                        
                        tmpFieldMean = [tmpFieldMean mean(pullBins)];
                        tmpSize = [tmpSize length(pfInds)*convFact];
                    end %pf
                    
                    peakFieldFirRate{g,1} = [peakFieldFirRate{g,1} mean(tmpPeak)];
                    meanFieldFirRate{g,1} = [meanFieldFirRate{g,1} mean(tmpFieldMean)];
                    
                    meanFirRate{g,1} = [meanFirRate{g,1} mean(uRateMaps(:))];
                    fieldSize{g,1} = [fieldSize{g,1} mean(tmpSize)];
                    
                end %if this unit has a place field
            end %unit
        end %day
    end %rat
end %group

cd(saveDir)
keyboard

%% SPATIAL CORRELATION

avgSpatCorr = cellfun(@mean, spatCorr);
semSpatCorr = cellfun(@semfunct, spatCorr);

figtitle = ['SpatialCorrelation_allCells'];

figure('Name', figtitle)

lh = NaN(1,2);
for g = 1:2
    hold on;
    lh(g) =  errorbar(1:4, avgSpatCorr(g,:), semSpatCorr(g,:), 'Color', cols{g}, 'LineWidth', 1); %just get the line handle for the first one for the legend
end %group

ylabel('Spatial correlation')
ylim([0 1])
xlim([0.5 4.5])
xticks(1:4)
xticklabels(combosTxt)
xlabel('Sessions compared')
legend(lh, groupNames, 'Location', 'southeast')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% RATE OVERLAP

avgRateOverlap = cellfun(@mean, rateOverlap);
semRateOverlap = cellfun(@semfunct, rateOverlap);

figtitle = ['RateOverlap_allCells'];

figure('Name', figtitle)

lh = NaN(1,2);
for g = 1:2
    hold on;
    lh(g) =  errorbar(1:4, avgRateOverlap(g,:), semRateOverlap(g,:), 'Color', cols{g}, 'LineWidth', 1); %just get the line handle for the first one for the legend
    
end %group

ylabel('Rate overlap')
ylim([0 1])
xlim([0.5 4.5])
xticks(1:4)
xticklabels(combosTxt)
xlabel('Sessions compared')
legend(lh, groupNames, 'Location', 'southeast')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% SPAT INFO

figtitle = 'SpatialInformation';
figure('Name', figtitle)

yLab = {'Spatial information (bits/spike)'};
bar_and_dotplot(spatInfo, yLab, groupNames, cols);

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% SPARSITY

figtitle = 'SpatialSparsity';
figure('Name', figtitle)

yLab = {'Spatial sparsity'};
bar_and_dotplot(sparsity, yLab, groupNames, cols);

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% PLACE FIELD SIZE

figtitle = 'PlaceFieldSize';
figure('Name', figtitle)

yLab = {'Place field size (cm)'};
bar_and_dotplot(fieldSize, yLab, groupNames, cols);

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% FIELD PER CELL

figtitle = 'PlaceFieldPerCell';
figure('Name', figtitle)

yLab = {'Place fields/cell'};
bar_and_dotplot(fieldPerCell, yLab, groupNames, cols);

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% MEAN FIRING RATE

figtitle = 'MeanFiringRate';
figure('Name', figtitle)

yLab = {'Mean firing rate (Hz)'};
bar_and_dotplot(meanFirRate, yLab, groupNames, cols);

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% PEAK FIRING RATE

figtitle = 'PeakPlaceFieldFiringRate';
figure('Name', figtitle)

yLab = {'Peak in-field firing rate (Hz)'};
bar_and_dotplot(peakFieldFirRate, yLab, groupNames, cols);

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% AVG IN-FIELD FIRING RATE

figtitle = 'AvgPlaceFieldFiringRate';
figure('Name', figtitle)

yLab = {'Average in-field firing rate (Hz)'};
bar_and_dotplot(meanFieldFirRate, yLab, groupNames, cols);

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% STATS?

if prepForStats == 1
    
    statSpatCorr = [];
    statRateOverlap = [];
    statSpatInfo = [];
    statFieldSize = [];
    statMeanFirRate = [];
    statPeakFieldFirRate = [];
    statMeanFieldFirRate = [];
    statFieldPerCell = [];
    
    for g = 1:2
        c = 1; %temporarily
        numU = length(spatCorr{g,c});
        for u = 1:numU
            tmpStatSC = g;
            tmpStatRO = g;
            
            for c = 1:4
                tmpStatSC = [tmpStatSC spatCorr{g,c}(u)];
                tmpStatRO = [tmpStatRO rateOverlap{g,c}(u)];
            end %combos
            statSpatCorr = [statSpatCorr; tmpStatSC];
            statRateOverlap = [statRateOverlap; tmpStatRO];
            
        end %unit
        
        for u = 1:length(spatInfo{g})
            statSpatInfo = [statSpatInfo; g spatInfo{g}(u)];
            
            statFieldSize = [statFieldSize; g fieldSize{g}(u)];
            
            statMeanFirRate = [statMeanFirRate; g meanFirRate{g}(u)];
            statPeakFieldFirRate = [statPeakFieldFirRate; g peakFieldFirRate{g}(u)];
            statMeanFieldFirRate = [statMeanFieldFirRate; g meanFieldFirRate{g}(u)];
            
            statFieldPerCell = [statFieldPerCell fieldPerCell{g}(u)];
        end %u - spat info
    end %group
    keyboard
end %prep for stats

end %function