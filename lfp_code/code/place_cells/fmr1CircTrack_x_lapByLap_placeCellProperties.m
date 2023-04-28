function fmr1CircTrack_x_lapByLap_placeCellProperties(group)
% function fmr1CircTrack_x_lapByLap_placeCellProperties(group)
% 
% PURPOSE:
%   Plot place cell propertes across laps to visualize stability across the
%   day.
% 
% INPUT:
%   group = data struct.
% 
% OUTPUT:
%   F1: Place cell firing rate across laps (in-field and peak).
%   F2: Place field center.
%   F3: Place field size ratio.
%   F4: Place field skewness. (NOTE: something is wrong with this).
%   F5: Place field firing rate asymmetry index (FRAI).
% 
% 
% MMD
% 10/2021
% Colgin Lab

%% OPTIONS
saveOrNot = 0;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\placeCellProperties\byLap';

prepForStats = 0;

minSpkPerLap = 3; %min # of in-field spikes per lap to be included - 3 in Feng Silva & Foster 2015, 1 in Mehta et al. 1997
minSpkBinsPerLap = 2; %as above

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
%% INITIALIZE

spatBinSz = 4;
velFilt = 1;
durCrit = 1;

lapFirRate = cell(2,numLaps);
bFirRate = cell(2,1);
lapPkFir = cell(2,numLaps);
bPkFir = cell(2,1);
lapPfSize = cell(2,numLaps);
lapPfCent = cell(2,numLaps);
lapPfSizeRatio = cell(2,numLaps);
lapPfSkew = cell(2,numLaps);
lapFRAI = cell(2,numLaps);

degCmConv = (pi*100)/360; %track has 1 m diameter

cols = {'Blue', 'Red'};

curDir = pwd;
cd(saveDir)

%% GET DATA


for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                fprintf('\t\t\tUnit %d/%d\n', u, length(group(g).rat(r).day(d).xAllBeginUnitInfo))
                tetNum = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID(1);
                clustNum = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID(2);
                
                xBegRateMap = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap;
                xBegPf = group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf;
                
                for p = 1:length(xBegPf)
                    passCntr = 0;
                    
                    pfRadPos = xBegPf(p).radPos;
                    startField = pfRadPos(1);
                    endField = pfRadPos(end);
                    pfLen = abs(rad2deg(circ_dist(deg2rad(startField), deg2rad(endField))));  %circular distance - pf can cross zero
                    pfLenCm = pfLen * degCmConv;
                    pfCent = wrapTo360(startField + pfLen/2);
                    
                    if endField > startField
                        crossZero = 0;
                    else
                        crossZero = 1;
                    end %determine if pf crosses 0
                    
                    rMapByPass = zeros(numLaps,length(xBegPf(p).inds));
                    pfPasses = cell(4,1); %initialize
                    
                    for b = 1:4
                        radPos = group(g).rat(r).day(d).begin(b).radPos;
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                        
                        pfPassBnry = zeros(1,size(radPos,1));
                        if crossZero == 0
                            pfPassBnry(radPos(:,2)>=startField & radPos(:,2)<=endField) = 1;
                        else
                            pfPassBnry(radPos(:,2)>=startField & radPos(:,2)<=360) = 1; %from start-360 and 0-end
                            pfPassBnry(radPos(:,2)>=0 & radPos(:,2)<=endField) = 1;
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
                            
                            if sum(spkCnts(xBegPf(p).inds)) < 3 %at least minSpkPerLap spikes
                                stableCheck = 0;
                                break %out of pass loop
                            end %at least 3 spikes
                            if length(find(spkCnts(xBegPf(p).inds))) < minSpkBinsPerLap %in at least minSpkBinsPerLap different position bins
                                stableCheck = 0;
                                break %out of pass loop
                            end %at least 2 pos bins
                            if passCntr <= numLaps %save the ratemap for the first 10 passes - for getting firing rate
                                rMapByPass(passCntr,:) = tmpMap(xBegPf(p).inds);
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
                        continue
                    end %continue to next unit
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
                        radPos = group(g).rat(r).day(d).begin(b).radPos;
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                        
                        for ps = 1:size(pfPasses{b},1)
                            passCntr = passCntr + 1;
                            if passCntr > numLaps
                                break %out of begin loop
                            end %more laps than we need for fig
                            cutSpkTms = spkTms(spkTms >= pfPasses{b}(ps,1) & spkTms <=pfPasses{b}(ps,2)); %pull out spike time just from this pass
                            psSpkBins = find(rMapByPass(passCntr,:));
                            psPfLen = abs(rad2deg(circ_dist(deg2rad(xBegPf(p).radPos(psSpkBins(1))), deg2rad(xBegPf(p).radPos(psSpkBins(end))))));
                            psPfCent = wrapTo360(xBegPf(p).radPos(psSpkBins(1)) + psPfLen/2);
                            tmpDist = rad2deg(circ_dist(deg2rad(pfCent), deg2rad(psPfCent)));
                            pfSkew = skewness(rMapByPass(passCntr,:))/std(rMapByPass(passCntr,:))^3;
                            
                            halfSpks = ceil(length(cutSpkTms)/2);
                            FR1 = halfSpks/(cutSpkTms(halfSpks) - cutSpkTms(1));
                            FR2 = halfSpks/(cutSpkTms(end) - cutSpkTms(halfSpks));
                            psFRAI = (FR1-FR2)/(FR1+FR2);
                            
                            lapFirRate{g,passCntr} = [lapFirRate{g,passCntr} mean(rMapByPass(passCntr,:))];
                            lapPkFir{g,passCntr} = [lapPkFir{g,passCntr} max(rMapByPass(passCntr,:))];
                            
                            
                            lapPfCent{g,passCntr} = [lapPfCent{g,passCntr} tmpDist*degCmConv]; %in cm
                            lapPfSizeRatio{g,passCntr} = [lapPfSizeRatio{g,passCntr} psPfLen/pfLen];
                            lapPfSkew{g,passCntr} =[lapPfSkew{g,passCntr} pfSkew];
                            lapFRAI{g,passCntr} = [lapFRAI{g,passCntr} psFRAI];
                            if isnan(psFRAI)
                                keyboard
                            end
                            
                            if prepForStats == 1
                                pfStatFr(passCntr+1) = mean(rMapByPass(passCntr,:));
                                pfStatPkFr(passCntr+1) = max(rMapByPass(passCntr,:));
                                pfStatCent(passCntr+1) = tmpDist*degCmConv;
                                pfStatSizeRat(passCntr+1) = psPfLen/pfLen;
                                pfStatSkew(passCntr+1) = pfSkew;
                                pfStatFRAI(passCntr+1) = psFRAI;
                            end %stats
                        end %passes in this begin
                    end %begin
                    bFirRate{g} = [bFirRate{g} mean(xBegRateMap(xBegPf(p).inds))];
                    bPkFir{g} = [bPkFir{g} max(xBegRateMap(xBegPf(p).inds))];
                end %pass
            end %unit
        end %day
    end %rat
end %group

keyboard

%% FIG 1 - PLACE CELL FIR RATE

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
ylim([0 6.5])
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


%% FIG 2 - PLACE FIELD CENTER

figtitle = 'PlaceFieldCenter_acrossLaps';

figure('Name', figtitle, 'Position', [594 467 699 420])

meanCent = cellfun(@mean, lapPfCent);
semCent = cellfun_emptyCells(@semfunct, lapPfCent, 0);

lh = nan(1,2);
leg = cell(2,1);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanCent(g,:), semCent(g,:), semCent(g,:), 'Color', rgb(cols{g}));
    
    leg{g} = [group(g).name ' n = ' num2str(length(lapPfCent{g})) ' cells'];
end %group

xlim([0 numLaps+1])
xlabel('Laps')
xticks([1:numLaps])
% xticklabels({[1:numLaps], 'Entire day'})
zero_line;

ylabel('Place field center (cm)')
% ylim([0 0.6])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 3 - PF SIZE RATIO

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

%% FIG 4 - SKEWNESS

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