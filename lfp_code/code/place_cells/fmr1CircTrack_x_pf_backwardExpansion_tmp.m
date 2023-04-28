function fmr1CircTrack_x_pf_backwardExpansion_tmp(group)
% function fmr1CircTrack_x_pf_backwardExpand(group)
%
% PURPOSE:
%   The purpose of this function is to look at the backwards expansion of
%   place cells in WT and KO rats.


%% OPTIONS

% wholeDay = 1; %to use the pf from the whole day (1) to just across the first min number of laps
minSpkPerLap = 3; %min # of in-field spikes per lap to be included - 3 in Feng Silva & Foster 2015, 1 in Mehta et al. 1997
minSpkBinsPerLap = 2; %as above

prepForStats = 0;


%% INITIALIZE
degCmConv = (100*pi)/360;
degBinCtrs = 2:4:360;

numBins = 100; %for normalizing

runThresh = 5; %cm/s

spMap = [1:2; 3:4; 5:6];
histMap = [5 6];
bMap = [1 4];
lpNames = {'first', 'last'};
tmpCols = {'Gray', 'Black'};

spatBinSz = 4; %degrees
velFilt = 1;
durCrit = 1;

frxPf = cell(2,2); %for fig 1

minLapDay = inf;
for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            dLaps = 0;
            for b = 1:4
                dLaps = dLaps + size(group(g).rat(r).day(d).begin(b).lapTms,1);
            end %begin
            if dLaps < minLapDay
                minLapDay = dLaps;
            end
        end %day
    end %rat
end %group
numLaps = minLapDay - 1; %number of laps to plot

lapFirRate = cell(2,numLaps);
bFirRate = cell(2,1);
lapPkFir = cell(2,numLaps);
bPkFir = cell(2,1);
lapPfSize = cell(2,numLaps);
lapPfCent = cell(2,numLaps);
bPfCent = cell(2,1);
lapPfSizeRatio = cell(2,numLaps);
lapPfSkew = cell(2,numLaps);
lapFRAI = cell(2,numLaps);

cols = {'Blue', 'Red'};

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
            %             for lpInd = 1:2
            %                 begRadPos{lpInd} = group(g).rat(r).day(d).begin(bMap(lpInd)).radPos;
            %                 begCoords{lpInd} = group(g).rat(r).day(d).begin(bMap(lpInd)).coords;
            %                 begSpd{lpInd} = smooth_runspeed(get_runspeed(begCoords{lpInd}));
            %             end %lp ind
            
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                xBegPf = group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf;
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
                    
                    if crossZero == 0
                        pfRm = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap(xBegPf(p).inds(1):xBegPf(p).inds(end));
                    else
                        pfRm = [group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap(1:end) group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap(1:pfInds(end))];
                    end %whether cross 0
                    
                    pfLen = rad2deg(circ_dist(deg2rad(pfPos(end)), deg2rad(pfPos(1))));
                    
                    props = regionprops(true(size(pfRm)),  pfRm, 'WeightedCentroid'); %get center of mass
                    pfCentNorm = props.WeightedCentroid(1) / length(pfRm); %normalize by the size of the pf
                    pfCent = wrapTo360(pfPos(1) + pfLen*pfCentNorm);
%                     pfHalfDist = abs(rad2deg(circ_dist(deg2rad(xBegPf(p).radPos(1)), deg2rad(xBegPf(p).radPos(end))))) / 2;

                    
                    rMapByPassPf = zeros(numLaps,length(xBegPf(p).inds));
                    rMapByPass = zeros(numLaps,length(degBinCtrs));
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
                            [tmpMap,~,~,spkCnts] = get_ratemap_circtrack(psSpkTms, psCoords, psRadPos, spatBinSz, velFilt, durCrit);
                            tmpMap = smooth_circtrack_ratemap(tmpMap, spatBinSz);
                            
                            if sum(spkCnts(xBegPf(p).inds)) < 3 %at least minSpkPerLap spikes
                                stableCheck = 0;
                                break %out of pass loop
                            end %at least 3 spikes
                            if length(find(spkCnts(xBegPf(p).inds))) < minSpkBinsPerLap %in at least minSpkBinsPerLap different position bins
                                stableCheck = 0;
                                break %out of pass loop
                            end %at least 2 pos bins
                            if passCntr <= numLaps %save the ratemap for the first 10 passes - for getting firing rate
                                rMapByPassPf(passCntr,:) = tmpMap(xBegPf(p).inds);
                                rMapByPass(passCntr,:) = tmpMap;
                            end %include this pass
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
                    
                    %start with getting over all data ("b...")
                    bFirRate{g} = [bFirRate{g} mean(pfRm)];
                    bPkFir{g} = [bPkFir{g} xBegPf(p).pkFr];
                    bPfCent{g} = [bPfCent{g} pfCent];
                    
                    if prepForStats == 1
                        pfStatSlope = [g nan(1,numLaps)]; %initialize
                        pfCmStatSlope = [g nan(1,numLaps)];
                        pfStatPO = [g nan(1,numLaps)];
                        pfStatPR = [g nan(1,numLaps)];
                        pfStatAltPR = [g nan(1,numLaps)];
                        pfStatR2 = [g nan(1,numLaps)];
                        
                        pfStatFr = [g nan(1,numLaps)];
                        pfStatPkFr = [g nan(1,numLaps)];
                        pfStatCent = [g nan(1,numLaps)];
                        pfStatSizeRat = [g nan(1,numLaps)];
                        pfStatSkew = [g nan(1,numLaps)];
                        pfStatFRAI = [g nan(1,numLaps)];
                    end %stats
                    
                    passCntr = 0; %re-initialize
                    for b = 1:4
                        radPos = group(g).rat(r).day(d).begin(b).radPos; %for entire begin
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                        
                        for ps = 1:size(pfPasses{b},1)
                            passCntr = passCntr + 1;
                            if passCntr > numLaps
                                break %out of begin loop
                            end %more laps than we need for fig
                       minPkFr = 0.25;
                       minPfLen = 5;
                            psPf = get_circtrack_pfs(rMapByPass(passCntr,:), spatBinSz, minPkFr, minPfLen)
                            if length(psPf) > 1
                                keyboard
                            end
                            for pp = 1:length(psPf)
                                checkSame = intersect(xBegPf(p).inds, psPf(pp).inds);
                                if length(checkSame) < 5
                                    keyboard
                                end
                            end %pf in pass
                            try
                             checkDiff = diff(psPf(pp).inds);
                            catch; keyboard; end
                            crossZero = 0; %initialize as doesn't
                            if ~isempty(find(checkDiff~=1))
                                crossZero = 1;
                            end %check
                            
                            if crossZero == 0
                                psPfRm = rMapByPass(passCntr,psPf(pp).inds(1):psPf(pp).inds(end));
                            else
                                keyboard
%                                 psPfRm = [rMapByPass(passCntr,:)(psPf(pp).inds(1):end) rMapByPass(passCntr,:)(1:psPf(pp))];
                            end %whether cross 0
                            
                            psPfLen = rad2deg(circ_dist(deg2rad(psPf(pp).radPos(end)), deg2rad(psPf(pp).radPos(1))));
                            
                            props = regionprops(true(size(psPfRm)),  psPfRm, 'WeightedCentroid'); %get center of mass
                            pfCentNorm = props.WeightedCentroid(1) / length(psPfRm); %normalize by the size of the pf
                            psPfCent = wrapTo360(psPf(pp).radPos(1) + psPfLen*pfCentNorm);
                            centDiff = rad2deg(circ_dist(deg2rad(psPfCent), deg2rad(pfCent))); %get ditance of the current center to other pf center
                            
                             psSkew = skewness(psPfRm);
                            
                            psSpkTms = spkTms(spkTms >= pfPasses{b}(ps,1) & spkTms <=pfPasses{b}(ps,2)); %pull out spike time just from this pass of the overall field
                            psFRAI = get_frai(psSpkTms);
                            
                            lapFirRate{g,passCntr} = [lapFirRate{g,passCntr} mean(psPfRm)];
                            lapPkFir{g,passCntr} = [lapPkFir{g,passCntr} psPf(pp).pkFr];
                            
                            lapPfCent{g,passCntr} = [lapPfCent{g,passCntr} centDiff];
                            lapPfSizeRatio{g,passCntr} = [lapPfSizeRatio{g,passCntr} psPfLen/pfLen];
                            lapPfSkew{g,passCntr} =[lapPfSkew{g,passCntr} psSkew];
                            lapFRAI{g,passCntr} = [lapFRAI{g,passCntr} psFRAI];
                            if isnan(psFRAI)
                                keyboard
                            end

                            if prepForStats == 1
                                pfStatFr(passCntr+1) = mean(rMapByPassPf(passCntr,:));
                                pfStatPkFr(passCntr+1) = max(rMapByPassPf(passCntr,:));
                                pfStatCent(passCntr+1) = tmpDist*degCmConv;
                                pfStatSizeRat(passCntr+1) = psPfLen/pfLen;
                                pfStatSkew(passCntr+1) = psSkew;
                                pfStatFRAI(passCntr+1) = psFRAI;
                            end %stats
                        end %passes in this begin
                    end %begin
                    
                    for lpInd = 1:2
                        subplot(3,2,spMap(lpInd,g))
                        
                        b = bMap(lpInd);
                        if lpInd == 1
                            lp = 1;
                        else
                            lp = size(group(g).rat(r).day(d).begin(b).lapTms,1);
                        end %which lap in this begin
                        lpTms = group(g).rat(r).day(d).begin(b).lapTms(lp,:);
                        radPos = group(g).rat(r).day(d).begin(b).radPos;
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        runSpd = smooth_runspeed(get_runspeed(coords));
                        
                        lpPos = radPos(match(lpTms(1), radPos(:,1)):match(lpTms(2), radPos(:,1)),:);
                        lpCoords = coords(match(lpTms(1), coords(:,1)):match(lpTms(2), coords(:,1)),:);
                        spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms; %all spkTms from this begin
                        psSpkTms = spkTms(spkTms >= lpTms(1) & spkTms <= lpTms(2)); %spikes during the first lap
                        lpRateMap = get_ratemap_circtrack(psSpkTms, lpCoords, lpPos, spatBinSz, velFilt, durCrit);
                        lpRateMap = smooth_circtrack_ratemap(lpRateMap, spatBinSz);
                        if crossZero == 0
                            pfCut = lpRateMap(pfInds(1):pfInds(end));
                        else
                            pfCut = [lpRateMap(pfInds(1):end) lpRateMap(1:pfInds(end))];
                        end %whether or not cross zero
                        tmp = round(rescale(1:length(pfInds))* 100);
                        
                        resizeMap = zeros(1,numBins);
                        startBin = 1;
                        for bin = 1:length(tmp)-1
                            endBin = tmp(bin+1);
                            try
                                resizeMap(startBin:endBin) = pfCut(bin);
                            catch; keyboard; end
                            startBin = tmp(bin+1)+1;
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
                                        line([normPos  normPos], [yVal-.4 yVal+.4], 'Color', rgb(tmpCols{lpInd}))
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
coM = zeros(2,2);
for g = 1:2
    lh = zeros(1,2);
    for lpInd = 1:2
        subplot(3,2,spMap(lpInd,g))
        ylim([1 size(frxPf{g,lpInd},1)])
        %           line([0.5  0.5], [0 size(frxPf{g,lpInd},1)], 'Color', rgb(tmpCols{lpInd}), 'LineStyle', '--')
        
        subplot(3,2,histMap(g))
        hold on;
        lh(lpInd) = plot(linspace(0, 1, numBins), meanFrxPf{g,lpInd}, 'Color', rgb(tmpCols{lpInd}));
        
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
          line([coM(g,lpInd) coM(g,lpInd)], [0 size(frxPf{g,lpInd},1)], 'Color', rgb(tmpCols{lpInd}), 'LineStyle', '--')
    end %lap ind
    
    subplot(3,2,histMap(g))
    ylim([0 yMax])
    for lpInd = 1:2
        line([coM(g,lpInd) coM(g,lpInd)], [0 yMax], 'Color', rgb(tmpCols{lpInd}), 'LineStyle', '--')
    end %lap ind
end %group

%% FIG 2 - FIRING RATE PLOTS

%avg in-field firing rate
figtitle = 'PlaceCellFirRate_acrossLaps';

figure('Name', figtitle, 'Position', [248 486 1408 420])

subplot(1,2,1)
meanFR = cellfun(@mean, lapFirRate);
bMeanFr = cellfun(@mean, bFirRate);
semFR = cellfun(@semfunct, lapFirRate);
bsemFr = cellfun(@semfunct, bFirRate);

lh = nan(1,2);
leg = cell(2,1);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanFR(g,:), semFR(g,:), semFR(g,:), 'Color', rgb(cols{g}));
    errorbar(numLaps+2, bMeanFr(g), bsemFr(g), bsemFr(g), 'Color', rgb(cols{g}))
    leg{g} = [group(g).name ' n = ' num2str(length(lapFirRate{g})) ' cells'];
end %group

xlim([0 numLaps+3])
xlabel('Laps')
xticks([1:numLaps numLaps+2])
xticklabels({[1:numLaps], 'Entire day'})

ylabel('In-field firing rate (Hz)')
ylim([0 6])
legend(lh, leg, 'Location', 'northeastoutside')

subplot(1,2,2)
meanFR= cellfun(@mean, lapPkFir);
bMeanFr = cellfun(@mean, lapPkFir);
semFR = cellfun(@semfunct, lapPkFir);
bsemFr = cellfun(@semfunct, lapPkFir);

lh = nan(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanFR(g,:), semFR(g,:), semFR(g,:), 'Color', rgb(cols{g}));
    errorbar(numLaps+2, bMeanFr(g), bsemFr(g), bsemFr(g), 'Color', rgb(cols{g}))
end %group

xlim([0 numLaps+3])
xlabel('Laps')
xticks([1:numLaps numLaps+2])
xticklabels({[1:numLaps], 'Entire day'})

ylabel('Peak firing rate (Hz)')
ylim([0 23])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 3 - PLACE FIELD CENTER

figtitle = 'PlaceFieldCenter_acrossLaps';
figure('Name', figtitle, 'Position', [594 467 699 420])

meanCent = cellfun(@mean, lapPfCent);
bMeanCent = cellfun(@mean, bPfCent);
semCent = cellfun(@semfunct, lapPfCent);
bSemCent = cellfun(@semfunct, bPfCent);

lh = nan(1,2);
leg = cell(2,1);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanCent(g,:), semCent(g,:), semCent(g,:), 'Color', rgb(cols{g}));
      errorbar(numLaps+2, bMeanCent(g), bSemCent(g), bSemCent(g), 'Color', rgb(cols{g}))
    leg{g} = [group(g).name ' n = ' num2str(length(lapPfCent{g})) ' cells'];
end %group

xlim([0 numLaps+3])
xlabel('Laps')
xticks([1:numLaps numLaps+2])
xticklabels({[1:numLaps], 'Entire day'})

ylabel('Place field center (normalized distance in pf)')
% ylim([0 0.6])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 4 - PF SIZE RATIO

figtitle = 'PlaceFieldSizeRatio_acrossLaps';
figure('Name', figtitle, 'Position', [594 467 699 420])

meaRatio = cellfun(@mean, lapPfSizeRatio);
semRatio = cellfun(@semfunct, lapPfSizeRatio);

lh = nan(1,2);
leg = cell(2,1);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meaRatio(g,:), semRatio(g,:), semRatio(g,:), 'Color', rgb(cols{g}));
    
    leg{g} = [group(g).name ' n = ' num2str(length(lapPfSizeRatio{g})) ' cells'];
end %group

xlim([0 numLaps+1])
xlabel('Laps')
xticks([1:numLaps])

ylabel('Place field size ratio')
ylim([0.5 0.8])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 5 - SKEWNESS

figtitle = 'PlaceFieldSkewness_acrossLaps';
figure('Name', figtitle, 'Position', [594 467 699 420])

meanSkew= cellfun(@mean, lapPfSkew);
semSkew = cellfun(@semfunct, lapPfSkew);

lh = nan(1,2);
leg = cell(2,1);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanSkew(g,:), semSkew(g,:), semSkew(g,:), 'Color', rgb(cols{g}));
    
    leg{g} = [group(g).name ' n = ' num2str(length(lapPfSkew{g})) ' cells'];
end %group

xlim([0 numLaps+1])
xlabel('Laps')
xticks([1:numLaps])

ylabel('Skewness')

legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 5 - FRAI

figtitle = 'PlaceFieldFRAI_acrossLaps';
figure('Name', figtitle, 'Position', [594 467 699 420])

meanFRAI = cellfun(@mean, lapFRAI);
semFRAI = cellfun(@semfunct, lapFRAI);

lh = nan(1,2);
leg = cell(2,1);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanFRAI(g,:), semFRAI(g,:), semFRAI(g,:), 'Color', rgb(cols{g}));
    
    leg{g} = [group(g).name ' n = ' num2str(length(lapFRAI{g})) ' cells'];
end %group

xlim([0 numLaps+1])
xlabel('Laps')
xticks([1:numLaps])

ylabel('Firing rate aymmsetry index')
zero_line;
ylim([-0.4 0.2])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

cd(curDir)

end %function