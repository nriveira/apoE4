function adprojCircTrack_x_pf_backwardExpansion_firstLastLap(group)
% function fmr1CircTrack_x_pf_backwardExpand(group)
%
% PURPOSE:
%   The purpose of this function is to look at the backwards expansion of
%   place cells in WT and KO rats. This uses a comparison between the first
%   and last lap on the track. The only cells included are those that have
%   stable place fields (see code for definition).
% 
% INPUT:
%   group = data struct
% 
% OUTPUT:
%   F1: (Top) Raster plots of spikes by normalized distance in place field
%       for first and last laps for both groups. Dashed lines: center of
%       mass of mean firing rate distribution. (Bottom) Mean firing rate
%       distribution across normalized distance in place field for both
%       groups. Dashed lines: center of mass of mean firing rate
%       distribution.
%   F2: Place field center distance (degrees from center of mass from
%       average place field) for first and last place field pass for both
%       groups.
%   F3: Place field peak firing rate position distance (degrees from
%       peak firing position in average place field) for first and last
%       place field passes for both groups.
%   F4: Place field skewness for first and last place field pass for both
%       groups.
%   F5: Place field FRAI (firing rate asymmetry index) for first and last
%       place field pass for both groups.
% 
% MMD
% 2/2022
% Colgin Lab

%% OPTIONS

saveOrNot = 1;

minSpkPerLap = 3; %min # of in-field spikes per lap to be included - 3 in Feng Silva & Foster 2015, 1 in Mehta et al. 1997
minSpkBinsPerLap = 2; %as above

prepForStats = 0;

%% INITIALIZE

saveDir = 'E:\resultsFeb2023_AD_WT\backwardExpansion';
degCmConv = (100*pi)/360;
degBinCtrs = 2:4:360;

numBins = 100; %for normalizing

runThresh = 5; %cm/s

spMap = [1:2; 3:4; 5:6];
histMap = [5 6];
bMap = [1 4];
lpNames = {'first', 'last'};
tmpCols = {'LightSkyBlue', 'DarkBlue'; 'LightSalmon', 'FireBrick'};

spatBinSz = 4; %degrees
velFilt = 1;
durCrit = 1;

frxPf = cell(1,2); %for fig 1
centDiff = cell(1,2);
pkDiff = cell(1,2);
pfSkew = cell(1,2);
pfFRAI = cell(1,2);

curDir = pwd;
cd(saveDir)

%% START FIG 1

figtitle = 'PlaceField_spikeLocations';
figure('Name', figtitle, 'Position', [565 109 787 869])

for sp = 1:4
    subplot(3,2,sp)
    xlabel('Normalized distance in place field')
    ylabel('Cell number')
    if sp == 1
        title({group(1).name, 'First lap'})
    elseif sp == 2
        title({group(2).name, 'First lap'})
    else
        title('Last lap')
    end %which title
end %subplot initializing
for sp = 5:6
    subplot(3,2,sp)
    xlabel('Normalized distance in place field')
    ylabel('Firing rate (Hz)')
end %subplot initializing

%% GET DATA/FILL IN FIG 1

for g = 1:2
    yVal = 1; %initialize for this group
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
%                         for lpInd = 1:2
%                             begRadPos{lpInd} = group(g).rat(r).day(d).begin(bMap(lpInd)).radPos;
%                             begCoords{lpInd} = group(g).rat(r).day(d).begin(bMap(lpInd)).coords;
%                             begSpd{lpInd} = smooth_runspeed(get_runspeed(begCoords{lpInd}));
%                         end %lp ind
            
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
%                 xBegPf = group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf;
                xBegPf = get_circtrack_pfs_m2(group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap, group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap); %use method 2 code so as not to cut off the bakward expansion
                if isempty(xBegPf)
                    continue %to next unit
                end %no pf
                for p = 1:length(xBegPf)
                    checkDiff = diff(xBegPf(p).inds);
                    crossZero = 0; %initialize as doesn't
                    if ~isempty(find(checkDiff~=1))
                        crossZero = 1;
                    end %check
                    pfInds = xBegPf(p).inds;
                    pfPos = xBegPf(p).radPos;
                    
                    pfRm = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap(pfInds);
                    
%                     pfLen = rad2deg(circ_dist(deg2rad(pfPos(end)), deg2rad(pfPos(1))));
                    pfLen = (length(pfInds)-1)*spatBinSz;
                    if pfLen >= 180
                        continue
                    end
                    pfPkPos = xBegPf(p).pkPos;
                    
                    props = regionprops(true(size(pfRm)),  pfRm, 'WeightedCentroid'); %get center of mass
                    pfCentNorm = props.WeightedCentroid(1) / length(pfRm); %normalize by the size of the pf
                    pfCent = wrapTo360(pfPos(1) + pfLen*pfCentNorm);
                    
                    %                     rMapByPass = zeros(numLaps,length(xBegPf(p).inds));
                    %                     smrMapByPass = zeros(numLaps,length(xBegPf(p).inds));
                    pfPasses = cell(4,1); %initialize
                    passCntr = 0;
                    %time to check whether pf is stable across all passes
                    for b = 1:4
                        radPos = group(g).rat(r).day(d).begin(b).radPos;
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                        
                        pfPassBnry = zeros(1,size(radPos,1));
                        if crossZero == 0
                            pfPassBnry(radPos(:,2)>=pfPos(1) & radPos(:,2)<=pfPos(end)) = 1;
                        else
                            pfPassBnry(radPos(:,2)>=pfPos(1) & radPos(:,2)<=360) = 1; %from start-360 and 0-end
                            pfPassBnry(radPos(:,2)>=0 & radPos(:,2)<= pfPos(end)) = 1;
                        end %cross zero
                        
                       pfPassChunks = bwconncomp(pfPassBnry);
                        for c = 1:length(pfPassChunks.PixelIdxList)
                            tmpInds = pfPassChunks.PixelIdxList{c};
                            passDist = abs(rad2deg(circ_dist(deg2rad(radPos(tmpInds(1),2)), deg2rad(radPos(tmpInds(end),2)))));
                            % Make sure rat traverses almost the whole field
                            if passDist >= pfLen - 5 % a little room for bin rounding error
                                pfPasses{b} = [pfPasses{b}; radPos(tmpInds(1),1) radPos(tmpInds(end),1)];
                            end %through whole field
                        end %chunks
                        
                        stableCheck = []; %empty
                        
                        for ps = 1:size(pfPasses{b},1) %stable = fire minSpkPerLap (3?) spike in minSpkBinsPerLap (2?) dif bins on each lap
                            passCntr = passCntr + 1;
                            psSpkTms = spkTms(spkTms>= pfPasses{b}(ps,1) & spkTms<=pfPasses{b}(ps,2));
                            psRadPos = radPos(radPos(:,1)>=pfPasses{b}(ps,1) & radPos(:,1)<=pfPasses{b}(ps,2),:);
                            psCoords = coords(coords(:,1)>=pfPasses{b}(ps,1) & coords(:,1)<=pfPasses{b}(ps,2),:);
                            try
                            [tmpMap,~,~,spkCnts] = get_ratemap_circtrack(psSpkTms, psCoords, psRadPos, spatBinSz, velFilt, durCrit);
                            catch; keyboard; end
                            %                             tmpSm = smooth_circtrack_ratemap(tmpMap, spatBinSz);
                            if sum(spkCnts(xBegPf(p).inds)) < minSpkPerLap %at least minSpkPerLap spikes
                                stableCheck = 0;
                                break %out of pass loop
                            end %at least 3 spikes
                            if length(find(spkCnts(xBegPf(p).inds))) < minSpkBinsPerLap %in at least minSpkBinsPerLap different position bins
                                stableCheck = 0;
                                break %out of pass loop
                            end %at least 2 pos bins
                            %                             if passCntr <= numLaps %save the ratemap for the first 10 passes - for getting firing rate
                            %                                 rMapByPass(passCntr,:) = tmpMap(xBegPf(p).inds);
                            %                                 smrMapByPass(passCntr,:) = tmpSm(xBegPf(p).inds);
                            %                             end %include this pass
                            if b == 4 && ps == size(pfPasses{b},1) %if this is the last pass and we haven't broken yet
                                stableCheck = 1;
                            end %got through all passes
                        end %pass
                        
                        if stableCheck == 0
                            break %out of begin loop
                        end %stable check
                    end %begin
                    
                    if stableCheck == 0
                        continue %to next pf
                    end %not stable
                    
                    %okay! if we have reached this point, this pf is "stable". we can move onto getting the data
                    %now to check for run speed on the first and last pass
                    spdCheck = 1; %initialize as good
                    for lpInd = 1:2
                        b = bMap(lpInd);
                        if lpInd == 1
                            ps = 1; %do it as passes instead of laps - first and last time rat passed through pf
                        else
                            ps = size(pfPasses{b},1);
                        end %which lap in this begin
                        try
                         psTms = pfPasses{b}(ps,:);
                        catch; keyboard; end
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        psCoords = coords(match(psTms(1), coords(:,1)):match(psTms(2), coords(:,1)),:);
                        runSpd = smooth_runspeed(get_runspeed(psCoords));
                        if min(runSpd) < runThresh
                            spdCheck = 0;
                        end %fails
                    end %lapInd
                    if spdCheck == 0
                        continue
                    end %doesn't pass
                    
                    for lpInd = 1:2
                        subplot(3,2,spMap(lpInd,g))
                        
                        b = bMap(lpInd);
                        if lpInd == 1
                            ps = 1; %do it as passes instead of laps - first and last time rat passed through pf
                        else
                            ps = size(pfPasses{b},1);
                        end %which lap in this begin
                        
                        psTms = pfPasses{b}(ps,:);
                        radPos = group(g).rat(r).day(d).begin(b).radPos;
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        runSpd = smooth_runspeed(get_runspeed(coords));
                        spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms; %all spkTms from this begin
                        
                        psPos = radPos(match(psTms(1), radPos(:,1)):match(psTms(2), radPos(:,1)),:);
                        psCoords = coords(match(psTms(1), coords(:,1)):match(psTms(2), coords(:,1)),:);
                        psSpkTms = spkTms(spkTms >= psTms(1) & spkTms <= psTms(2)); %spikes during the first lap
                        if isempty(psSpkTms)
                            continue
                        end
                        psRateMap = get_ratemap_circtrack(psSpkTms, psCoords, psPos, spatBinSz, velFilt, durCrit);
                        psRateMap = smooth_circtrack_ratemap(psRateMap, spatBinSz);
                        if crossZero == 0
                            psRm = psRateMap(pfInds(1):pfInds(end));
                        else
                            psRm = [psRateMap(pfInds(1):end) psRateMap(1:pfInds(end))];
                        end %whether or not cross zero
                        
                        [~, pfPkInd] = max(psRm);
                        psPkPos = wrapTo360(pfPkInd*spatBinSz + pfPos(1));
                        tmpPkDiff = rad2deg(circ_dist(deg2rad(psPkPos), deg2rad(pfPkPos)));
                        pkDiff{g,lpInd} = [pkDiff{g,lpInd} tmpPkDiff];
                        
                        props = regionprops(true(size(psRm)),  psRm, 'WeightedCentroid'); %get center of mass
                        psCentNorm = props.WeightedCentroid(1) / length(psRm); %normalize by the size of the pf
                        psCent = wrapTo360(pfPos(1) + pfLen*psCentNorm);
                        tmpCentDiff = rad2deg(circ_dist(deg2rad(psCent), deg2rad(pfCent)));
                        %                         if isnan(tmpCentDiff)
                        %                             keyboard
                        %                         end
                        centDiff{g,lpInd} = [centDiff{g,lpInd} tmpCentDiff];
                        
                        tmpSkew = skewness(psRm);
                        pfSkew{g,lpInd} = [pfSkew{g,lpInd} tmpSkew];
                        
                        tmpFRAI = get_frai(psSpkTms);
                        pfFRAI{g,lpInd} = [pfFRAI{g,lpInd} tmpFRAI];
                        
                        newTimeBins = round(rescale(1:length(pfInds))* 100);
                        
                        resizeMap = zeros(1,numBins);
                        startBin = 1;
                        for bin = 1:length(newTimeBins)-1
                            endBin = newTimeBins(bin+1);
                            try
                                resizeMap(startBin:endBin) = psRm(bin);
                            catch; keyboard; end
                            startBin = newTimeBins(bin+1)+1;
                        end %bin
                        
                        frxPf{g,lpInd} = [frxPf{g,lpInd}; resizeMap];
                        
                        for st = 1:length(psSpkTms)
                            spkPosInd = match(psSpkTms(st), radPos(:,1));
                            if runSpd(spkPosInd,2) >= runThresh
                                spkPos = radPos(spkPosInd,2);
                                
                                if crossZero == 0
                                    if spkPos >= pfPos(1) && spkPos <= pfPos(end)
                                        distFromStart = rad2deg(circ_dist(deg2rad(spkPos), deg2rad(pfPos(1))));
                                        normPos = distFromStart / pfLen; %
                                        if normPos > 1
                                            keyboard
                                        end
                                        line([normPos  normPos], [yVal-.4 yVal+.4], 'Color', rgb(tmpCols{g,lpInd}))
                                    end %within pf spike
                                end %cross zero
                                
                            end %speed check
                        end %spk tms
                        
                    end %lpInd
                    yVal = yVal + 1;
                end %place field
            end %unit
        end %day
    end %rat
end %group

keyboard
%% FINISH FIG 1
meanFrxPf = cellfun(@mean, frxPf, 'UniformOutput', false);

yMax = 0; %for group fr plots
coM = zeros(1,2);
for g = 1:2
    lh = zeros(1,2);
    for lpInd = 1:2
        subplot(3,2,spMap(lpInd,g))
        ylim([1 size(frxPf{g,lpInd},1)])
        %           line([0.5  0.5], [0 size(frxPf{g,lpInd},1)], 'Color', rgb(tmpCols{lpInd}), 'LineStyle', '--')
        
        subplot(3,2,histMap(g))
        hold on;
        lh(lpInd) = plot(linspace(0, 1, numBins), meanFrxPf{g,lpInd}, 'Color', rgb(tmpCols{g,lpInd}));
        
        props = regionprops(true(size(meanFrxPf{g,lpInd})),  meanFrxPf{g,lpInd}, 'WeightedCentroid'); %get center of mass
        coM(g,lpInd) = props.WeightedCentroid(1) / numBins; %divided by num bins for our xscale from 0 to 1
        
        yRange = get(gca, 'YLim');
        if yRange(2) > yMax
            yMax = yRange(2);
        end %fill in max
    end %lpInd
end %group

tmpMap = [1:2; 3:4];
for g = 1:2
    for lpInd = 1:2
        subplot(3,2,tmpMap(lpInd,g))
        line([coM(g,lpInd) coM(g,lpInd)], [0 size(frxPf{g,lpInd},1)], 'Color', rgb(tmpCols{g,lpInd}), 'LineStyle', '--')
    end %lap ind
    
    subplot(3,2,histMap(g))
    ylim([0 yMax])
    for lpInd = 1:2
        line([coM(g,lpInd) coM(g,lpInd)], [0 yMax], 'Color', rgb(tmpCols{g,lpInd}), 'LineStyle', '--')
    end %lap ind
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

% For the below figures, WeightedProportion function is not found. 
% Check again. 
%% FIG 2

figtitle = 'PlaceFieldCenter_firstLastLap';
figure('Name', figtitle, 'Position', [478 523 906 420])

minVal = floor(min(cat(2, centDiff{:})));
maxVal = ceil(max(cat(2, centDiff{:})));
sigma = 3;

for g = 1:2
    subplot(1,2,g)
    lh = zeros(1,2);
    leg = cell(1,2);
    for lpInd = 1:2
        hold on;
        [wDist, xAx] = WeightedProportion(centDiff{g,lpInd}', minVal, maxVal, sigma);
       lh(lpInd) = plot(xAx, wDist, 'Color', rgb(tmpCols{g,lpInd}));
       leg{lpInd} = [lpNames{lpInd} ' lap'];
    end %lap ind
    
    title([group(g).name ' (n = ' num2str(length(centDiff{g,lpInd})) ' cells)'])
    ylabel('Proportion of cells')
    xlabel('Place field center difference (deg)')
       legend(lh, leg, 'AutoUpdate','off', 'Location', 'northwest')
end %group

yAx = same_axes;
for g = 1:2
      subplot(1,2,g)
    for lpInd = 1:2
        line([mean(centDiff{g,lpInd}) mean(centDiff{g,lpInd})], [0 yAx(2)], 'Color', rgb(tmpCols{g,lpInd}), 'LineStyle', '--')
    end %lpInd
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 3

figtitle = 'PlaceFieldPkPos_firstLastLap';
figure('Name', figtitle, 'Position', [478 523 906 420])

minVal = floor(min(cat(2, pkDiff{:})));
maxVal = ceil(max(cat(2, pkDiff{:})));
sigma = 3;

for g = 1:2
    subplot(1,2,g)
    lh = zeros(1,1);
    leg = cell(1,1);
    for lpInd = 1:2
        hold on;
        [wDist, xAx] = WeightedProportion(pkDiff{g,lpInd}', minVal, maxVal, sigma);
       lh(lpInd) = plot(xAx, wDist, 'Color', rgb(tmpCols{g,lpInd}));
       leg{lpInd} = [lpNames{lpInd} ' lap'];
    end %lap ind
    
    title([group(g).name ' (n = ' num2str(length(pkDiff{g,lpInd})) ' cells)'])
    ylabel('Proportion of cells')
    xlabel('Place field peak position difference (deg)')
       legend(lh, leg, 'AutoUpdate','off', 'Location', 'northwest')
end %group

yAx = same_axes;
for g = 1:2
      subplot(1,2,g)
    for lpInd = 1:2
        line([mean(pkDiff{g,lpInd}) mean(pkDiff{g,lpInd})], [0 yAx(2)], 'Color', rgb(tmpCols{g,lpInd}), 'LineStyle', '--')
    end %lpInd
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 4

figtitle = 'PlaceFieldSkew_firstLastLap';
figure('Name', figtitle, 'Position', [478 523 906 420])

minVal = floor(min(cat(2, pfSkew{:})));
maxVal = ceil(max(cat(2, pfSkew{:})));
sigma = 0.5;

for g = 1:2
    subplot(1,2,g)
     lh = zeros(1,1);
    leg = cell(1,1);
    for lpInd = 1:2
        hold on;
        [wDist, xAx] = WeightedProportion(pfSkew{g,lpInd}', minVal, maxVal, sigma);
        lh(lpInd) = plot(xAx, wDist, 'Color', rgb(tmpCols{g,lpInd}));
          leg{lpInd} = [lpNames{lpInd} ' lap'];
    end %lap ind
    
    title([group(g).name ' (n = ' num2str(length(pkDiff{g,lpInd})) ' cells)'])
    ylabel('Proportion of cells')
    xlabel('Place field skewness')
     legend(lh, leg, 'AutoUpdate','off', 'Location', 'northwest')
end %group

yAx = same_axes;
for g = 1:2
      subplot(1,2,g)
    for lpInd = 1:2
        line([mean(pfSkew{g,lpInd}) mean(pfSkew{g,lpInd})], [0 yAx(2)], 'Color', rgb(tmpCols{g,lpInd}), 'LineStyle', '--')
    end %lpInd
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 5

figtitle = 'PlaceFieldFRAI_firstLastLap';
figure('Name', figtitle, 'Position', [478 523 906 420])

minVal = floor(min(cat(2, pfFRAI{:})));
maxVal = ceil(max(cat(2, pfFRAI{:})));
sigma = 0.1;

for g = 1:2
    subplot(1,2,g)
     lh = zeros(1,1);
    leg = cell(1,1);
    for lpInd = 1:2
        hold on;
        [wDist, xAx] = WeightedProportion(pfFRAI{g,lpInd}', minVal, maxVal, sigma);
        lh(lpInd) = plot(xAx, wDist, 'Color', rgb(tmpCols{g,lpInd}));
           leg{lpInd} = [lpNames{lpInd} ' lap'];
    end %lap ind
    
   title([group(g).name ' (n = ' num2str(length(pkDiff{g,lpInd})) ' cells)'])
    ylabel('Proportion of cells')
    xlabel('Place field FRAI')
     legend(lh, leg, 'AutoUpdate','off', 'Location', 'northwest')
end %group

yAx = same_axes;
for g = 1:2
      subplot(1,2,g)
    for lpInd = 1:2
        line([mean(pfFRAI{g,lpInd}) mean(pfFRAI{g,lpInd})], [0 yAx(2)], 'Color', rgb(tmpCols{g,lpInd}), 'LineStyle', '--')
    end %lpInd
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

cd(curDir)

end %function