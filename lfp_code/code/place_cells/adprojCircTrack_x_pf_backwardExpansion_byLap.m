%function adprojCircTrack_x_pf_backwardExpansion_byLap(group)
% function fmr1CircTrack_x_pf_backwardExpand(group)
%
% WORK IN PROGRESS

% PURPOSE:
%   The purpose of this function is to look at the backwards expansion of
%   place cells in rats. This function looks at the single-lap
%   place field across the first (min number) of laps in a single day.
%
% INPUT:
%   group = data struct
%
% OUTPUT:
%   F1: In-field and peak in-field firing rate across laps and average
%       across day.
%   F2: Place field center difference (in degrees - distance between
%       overall and single-lap) across laps.
%   F3: Place field peak position difference (in degrees - distance between
%       overall and single-lap) across laps.
%   F4: Place field size ratio(lap pfsize/overall pf size) across laps.
%   F5: Place field skewness across laps.
%   F6: FRAI (firing rate asymmetry index) across laps.
%
% Jayanth
% 2/2023
%
% MMD
% 2/2022
% Colgin Lab

clc; keep group;

%% OPTIONS

saveOrNot = 1;
%
minSpkPerLap = 3; %min # of in-field spikes per lap to be included - 3 in Feng Silva & Foster 2015, 1 in Mehta et al. 1997
minSpkBinsPerLap = 2; %as above

numLaps = 10;

prepForStats = 0;

%% INITIALIZE

saveDir = 'E:\resultsFeb2023_AD_WT\backwardExpansion\byLap';
degCmConv = (100*pi)/360;
degBinCtrs = 2:4:360;

runThresh = 5; %cm/s

spatBinSz = 4; %degrees
velFilt = 1;
durCrit = 1;

% minLapDay = inf;
% for g = 1:2
%     for r = 1:length(group(g).rat)
%         for d = 1:length(group(g).rat(r).day)
%             dLaps = 0;
%             for b = 1:4
%                 dLaps = dLaps + size(group(g).rat(r).day(d).begin(b).lapTms,1);
%             end %begin
%             if dLaps < minLapDay
%                 minLapDay = dLaps;
%             end
%         end %day
%     end %rat
% end %group

% numLaps = minLapDay - 1; %number of laps to plot
% maxLapDay = -inf;

lapFirRate = cell(2,numLaps);
bFirRate = cell(2,1);
lapPkFir = cell(2,numLaps);
bPkFir = cell(2,1);
lapPkPos = cell(2,numLaps);
lapPfSize = cell(2,numLaps);
lapPfCent = cell(2,numLaps);
lapPfSizeRatio = cell(2,numLaps);
lapPfSkew = cell(2,numLaps);
lapFRAI = cell(2,numLaps);

cols = {'Blue', 'Red'};

% lapFirRate = cell(1,numLaps);
% bFirRate = cell(1,1);
% lapPkFir = cell(1,numLaps);
% bPkFir = cell(1,1);
% lapPkPos = cell(1,numLaps);
% lapPfSize = cell(1,numLaps);
% lapPfCent = cell(1,numLaps);
% lapPfSizeRatio = cell(1,numLaps);
% lapPfSkew = cell(1,numLaps);
% lapFRAI = cell(1,numLaps);

curDir = pwd;
cd(saveDir)

%% GET DATA

for g = 1:2
    
    for r = 1:length(group(g).rat)
        %for d = 1:length(group(g).rat(r).day)
%         
%         if g == 1
%             d = 5;
%         else
%             d = 3;
%         end

        for d = 3

          %for b = 1:4
          for b = 1
            
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                %                 xBegPf = group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf;
                %xBegPf = get_circtrack_pfs(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap, spatBinSz);
                %xBegPf = get_circtrack_pfs_m2(group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap, group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap); %use method 2 code so as not to cut off the backward expansion
                
                xBegPf = get_circtrack_pfs_backwardExpansion(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap); % this is considering the whole day's rate map of the unit
                %meg shared this function, so check this function carefully and see if it's getting right bounds for the extension
                
                numFields(u,g) = length(xBegPf);
                
                if isempty(xBegPf)
                    continue %to next unit
                end %no pf

%                 if length(xBegPf)>1
%                     continue %to next unit
%                 end % more than 1 pf per unit
                
                for p = 1:length(xBegPf) % number of place fields
                    
                    checkDiff = diff(xBegPf(p).inds);
                    crossZero = 0; %initialize as the place field does not cross 360 degrees
                    
                    if ~isempty(find(checkDiff~=1))
                        crossZero = 1; % place field crosses 360 degrees
                    end %check
                    
                    pfInds = xBegPf(p).inds; % spatial bin numbers for the expanded place field
                    pfPos = xBegPf(p).radPos; % radial position in degrees for the extended place field
                    
                    pfRm = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap(pfInds);
                    
                    pfLen = abs(rad2deg(circ_dist(deg2rad(pfPos(end)), deg2rad(pfPos(1))))); % place field size in degrees
                    %pfLen = (length(pfInds)-1)*spatBinSz;
                    
                    % to check if the extended place field size is more
                    % than half of the track, if so do not consider them
                    if pfLen > 180
                        continue
                    end
                    
                    pfPkPos = xBegPf(p).pkPos;
                    
                    %props = regionprops(true(size(pfRm)),  pfRm, 'WeightedCentroid'); %get center of mass for all the laps combined during the day
                    tmpfRm = [];
                    tmpfRm = pfRm./sum(pfRm);
                    wmean = 0;
                    for i = 1:length(tmpfRm)
                        wmean = wmean + (i * tmpfRm(i));
                    end
                    pfcom = round(wmean); %get center of mass for all the laps combined during the day
                    pfcomBin = xBegPf(p).radPos(pfcom);

                    %pfCentNorm = props.WeightedCentroid(1) / length(pfRm); %normalize by the size of the pf

                    pfCentNorm = pfcom / length(pfRm); 
                    pfCent = wrapTo360(pfPos(1) + pfLen*pfCentNorm);
                    
                    %                     rMapByPass = zeros(numLaps,length(xBegPf(p).inds));
                    %                     smrMapByPass = zeros(numLaps,length(xBegPf(p).inds));
                    pfPasses = cell(4,1); %initialize
                    pfPassesInds = cell(4,1); % initialize
                    %passCntr = 0; %pass/lap counter
                    
                    %time to check whether pf is stable across all passes
                    %for b = 1:4
                    for b = 1
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
                                pfPassesInds{b} = [pfPassesInds{b}; tmpInds(1) tmpInds(end)];
                                pfPasses{b} = [pfPasses{b}; radPos(tmpInds(1),1) radPos(tmpInds(end),1)];
                            end %through whole field
                        end %chunks
                        
                        %                                                 stableCheck = []; %empty
                        %
                        %                                                 for ps = 1:size(pfPasses{b},1) %stable = fire minSpkPerLap (3?) spike in minSpkBinsPerLap (2?) dif bins on each lap
                        %                                                     passCntr = passCntr + 1;
                        %                                                     psSpkTms = spkTms(spkTms>= pfPasses{b}(ps,1) & spkTms<=pfPasses{b}(ps,2));
                        %                                                     psRadPos = radPos(radPos(:,1)>=pfPasses{b}(ps,1) & radPos(:,1)<=pfPasses{b}(ps,2),:);
                        %                                                     psCoords = coords(coords(:,1)>=pfPasses{b}(ps,1) & coords(:,1)<=pfPasses{b}(ps,2),:);
                        %                                                     [tmpMap,~,~,spkCnts] = get_ratemap_circtrack(psSpkTms, psCoords, psRadPos, spatBinSz, velFilt, durCrit);
                        %                                                     tmpSm = smooth_circtrack_ratemap(tmpMap, spatBinSz);
                        %                                                     if sum(spkCnts(xBegPf(p).inds)) < minSpkPerLap %at least minSpkPerLap spikes
                        %                                                         stableCheck = 0;
                        %                                                         break %out of pass loop
                        %                                                     end %at least 3 spikes
                        %                                                     if length(find(spkCnts(xBegPf(p).inds))) < minSpkBinsPerLap %in at least minSpkBinsPerLap different position bins
                        %                                                         stableCheck = 0;
                        %                                                         break %out of pass loop
                        %                                                     end %at least 2 pos bins
                        % %                                                     if passCntr <= numLaps %save the ratemap for the first 10 passes - for getting firing rate
                        % %                                                         rMapByPass(passCntr,:) = tmpMap(xBegPf(p).inds);
                        % %                                                         smrMapByPass(passCntr,:) = tmpSm(xBegPf(p).inds);
                        % %                                                     end %include this pass
                        %                                                     if b == 4 && ps == size(pfPasses{b},1) %if this is the last pass and we haven't broken yet
                        %                                                         stableCheck = 1;
                        %                                                     end %got through all passes
                        %                                                 end %pass
                        %
                        %                                                 if stableCheck == 0
                        %                                                     break %out of begin loop
                        %                                                 end %stable check
                    end %begin
                    
                    %                                         if stableCheck == 0
                    %                                             continue %to next pf
                    %                                         end %not stable
                    
                    %start with getting over all data ("b...")
                    bFirRate{g} = [bFirRate{g} mean(pfRm)];
                    bPkFir{g} = [bPkFir{g} xBegPf(p).pkFr];
                    
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
                    %for b = 1:4
                    for b = 1   
                        radPos = group(g).rat(r).day(d).begin(b).radPos; %for entire begin
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                        
                        for ps = 1:size(pfPasses{b},1)
                            passCntr = passCntr + 1;
                            if passCntr > numLaps
                                break %out of begin loop
                            end %more laps than we need for fig
                           

                            psSpkTms = spkTms(spkTms >= pfPasses{b}(ps,1) & spkTms <=pfPasses{b}(ps,2)); %pull out spike time just from this pass
                            %                             if isempty(psSpkTms) || length(psSpkTms) < 3
                            %                                 continue %to next pass
                            %                             end %no spikes in pf on this lap
                            
                            psRadPos = radPos(radPos(:,1)>=pfPasses{b}(ps,1) & radPos(:,1)<=pfPasses{b}(ps,2),:);
                            psCoords = coords(coords(:,1)>=pfPasses{b}(ps,1) & coords(:,1)<=pfPasses{b}(ps,2),:);
                            try
                                tmpMap = get_ratemap_circtrack(psSpkTms, psCoords, psRadPos, spatBinSz, velFilt, durCrit); %get pass ratemap
                            catch 
                                %keyboard; 
                            end
                            psRawRm = tmpMap(xBegPf(p).inds);
                            tmpSm = smooth_circtrack_ratemap(tmpMap, spatBinSz); %smooth the ratemap
                            psRateMap = tmpSm(xBegPf(p).inds); %smoothened ratemap of the unit within the place field bounds
                            if length(psSpkTms) < 3 || sum(psRawRm) == 0 || isempty(find(psRawRm))
                                continue %to next pass
                            end %disqualifying criteria
                            
                            %props = regionprops(true(size(psRateMap)),  psRateMap, 'WeightedCentroid'); %get center of mass
                            
                            pstmpfRm = [];
                            pstmpfRm = psRateMap./sum(psRateMap);
                            
                            wmean = 0;
                            for i = 1:length(pstmpfRm)
                                wmean = wmean + (i * pstmpfRm(i));
                            end
                            psPfcom = round(wmean); %get center of mass for all the laps combined during the day
                            psPfcomBin = xBegPf(p).radPos(psPfcom)/4;
                            %centDiff = deg2rad(pfcomBin - psPfcomBin);

                            %pfCentNorm = props.WeightedCentroid(1) / length(pfRm); %normalize by the size of the pf

                            psPfCentNorm = psPfcom / length(pfRm); % normalize by the size of the overall pf
                            %psPfCentNorm = psPfcom / length(pstmpfRm); % normalize by the size of the lap/pass pf

                            psPfCent = wrapTo360(pfPos(1) + pfLen*psPfCentNorm);
                            
                            %psPfCentNorm = props.WeightedCentroid(1) /length(pfRm); %normalize by the size of the overall pf
                            %psPfCent = wrapTo360(pfPos(1) + pfLen*psPfCentNorm);
                            
                            centDiff = rad2deg(circ_dist(deg2rad(psPfCent), deg2rad(pfCent))); %get distance of the current pass pf center to the overall day or begin pf center
                            
                            psSpkBins = find(psRawRm); %find bins that had spikes within the pass pf
                            
                            try
                                psPfLen = abs(rad2deg(circ_dist(deg2rad(xBegPf(p).radPos(psSpkBins(1))), deg2rad(xBegPf(p).radPos(psSpkBins(end))))));
                            catch; keyboard; end
                            
                            [~, psPkPosInd] = max(psRateMap);
                            psPkPos = xBegPf(p).radPos(psPkPosInd);
                            psPkPosNorm = psPkPosInd/length(psRateMap);
                            pkDiff = rad2deg(circ_dist(deg2rad(psPkPos), deg2rad(pfPkPos))); %get distance of the current pk to other pf pk
                            
                            pfSkew = skewness(psRateMap);
                            
                            psFRAI = get_frai(psSpkTms);
                            
                            lapFirRate{g,passCntr} = [lapFirRate{g,passCntr} mean(psRateMap)];
                            lapPkFir{g,passCntr} = [lapPkFir{g,passCntr} max(psRateMap)];
                            lapPkPos{g,passCntr} = [lapPkPos{g,passCntr} pkDiff];
                            %lapPkPos{g,passCntr} = [lapPkPos{g,passCntr} psPkPosNorm];
                            lapPfCent{g,passCntr} = [lapPfCent{g,passCntr} centDiff];
                            %lapPfCent{g,passCntr} = [lapPfCent{g,passCntr} psPfCentNorm];
                            lapPfSizeRatio{g,passCntr} = [lapPfSizeRatio{g,passCntr} psPfLen/pfLen];
                            lapPfSkew{g,passCntr} =[lapPfSkew{g,passCntr} pfSkew];
                            lapFRAI{g,passCntr} = [lapFRAI{g,passCntr} psFRAI];
                            if isnan(psFRAI)
                                keyboard
                            end
                            
                            %                             if prepForStats == 1
                            %                                 pfStatFr(passCntr+1) = mean(smrMapByPass(passCntr,:));
                            %                                 pfStatPkFr(passCntr+1) = max(smrMapByPass(passCntr,:));
                            %                                 pfStatCent(passCntr+1) = tmpDist*degCmConv;
                            %                                 pfStatSizeRat(passCntr+1) = psPfLen/pfLen;
                            %                                 pfStatSkew(passCntr+1) = pfSkew;
                            %                                 pfStatFRAI(passCntr+1) = psFRAI;
                            %                             end %stats
                        end %passes in this begin
                    end %begin
                end %place field
            end %unit
          end %begin
        end %day
    end %rat
end %group
keyboard

%% FIG 1 - FIRING RATE PLOTS

%avg in-field firing rate
figtitle = 'begin1E_PlaceCellFirRate';
figure('Name', figtitle, 'Position', [248 486 1408 420])

subplot(1,2,1)
meanFR = cellfun(@mean, lapFirRate);
bMeanFr = cellfun(@mean, bFirRate);
semFR = cellfun(@semfunct, lapFirRate);
bsemFr = cellfun(@semfunct, bFirRate);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanFR(g,:), semFR(g,:), semFR(g,:), 'Color', rgb(cols{g}));
    errorbar(numLaps+2, bMeanFr(g), bsemFr(g), bsemFr(g), 'Color', rgb(cols{g}))
    leg{g} = [group(g).name ' n = ' num2str(length(lapFirRate{g})) ' fields (from the first lap)'];
end %group

xlim([0 numLaps+3])
xlabel('Laps')
xticks([1:5:numLaps numLaps+2])
xticklabels({[1:5:numLaps], 'Entire day'})

ylabel('In-field firing rate (Hz)')
ylim([0 6])
legend(lh, leg, 'Location', 'northeastoutside')

%pk in field firing rate
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
xticks([1:5:numLaps numLaps+2])
xticklabels({[1:5:numLaps], 'Entire day'})

ylabel('Peak firing rate (Hz)')
ylim([0 15])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 2 - PLACE FIELD CENTER

figtitle = 'begin1E_PlaceFieldCenter';
figure('Name', figtitle, 'Position', [594 467 699 420])

meanCent = cellfun_emptyCells(@mean, lapPfCent, NaN);
semCent = cellfun_emptyCells(@semfunct, lapPfCent, 0);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanCent(g,:), semCent(g,:), semCent(g,:), 'Color', rgb(cols{g}));
    %lh(g) = errorbar(numLaps+2, bMeanCent(g), bSemCent(g), bSemCent(g), 'Color', rgb(cols{g}))
    leg{g} = [group(g).name ' n = ' num2str(length(lapPfCent{g})) ' fields (from the first lap)'];
end %group

xlim([0 numLaps+1])
xlabel('Laps')
xticks(0:1:numLaps)
%yline(deg2rad(pfcomBin), '--');
%zero_line;
% xticks([1:numLaps numLaps+2])
% xticklabels({[1:numLaps], 'Entire day'})

ylabel('Place field center (deg)')
% ylim([0 0.6])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 3 - PK POS

figtitle = 'begin1E_PlaceFieldPkPos';
figure('Name', figtitle, 'Position', [594 467 699 420])

meanPk = cellfun_emptyCells(@mean, lapPkPos, NaN);
semPk = cellfun_emptyCells(@semfunct, lapPkPos, NaN);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanPk(g,:), semPk(g,:), semPk(g,:), 'Color', rgb(cols{g}));
    %       errorbar(numLaps+2, bMeanCent(g), bSemCent(g), bSemCent(g), 'Color', rgb(cols{g}))
    leg{g} = [group(g).name ' n = ' num2str(length(lapPkPos{g})) ' fields (from the first lap)'];
end %group

xlim([0 numLaps+1])
xlabel('Laps')
xticks(1:1:numLaps)
%zero_line;
% xticks([1:numLaps numLaps+2])
% xticklabels({[1:numLaps], 'Entire day'})

ylabel('Place field peak pos (deg)')
% ylim([0 0.6])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 4 - PF SIZE RATIO

figtitle = 'begin1E_PlaceFieldSizeRatio';
figure('Name', figtitle, 'Position', [594 467 699 420])

meaRatio = cellfun(@mean, lapPfSizeRatio);
semRatio = cellfun(@semfunct, lapPfSizeRatio);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meaRatio(g,:), semRatio(g,:), semRatio(g,:), 'Color', rgb(cols{g}));
    
    leg{g} = [group(g).name ' n = ' num2str(length(lapPfSizeRatio{g})) ' fields'];
end %group

xlim([0 numLaps+1])
xlabel('Laps')
xticks(1:1:numLaps)

ylabel('Place field size ratio')
%ylim([0 0.8])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 5 - SKEWNESS

figtitle = 'PlaceFieldSkewness';
figure('Name', figtitle, 'Position', [594 467 699 420])

meanSkew = cellfun(@mean, lapPfSkew);
semSkew = cellfun(@semfunct, lapPfSkew);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanSkew(g,:), semSkew(g,:), semSkew(g,:), 'Color', rgb(cols{g}));
    
    leg{g} = [group(g).name ' n = ' num2str(length(lapPfSkew{g})) ' fields'];
end %group

xlim([0 numLaps+1])
xlabel('Laps')
xticks(1:5:numLaps)

ylabel('Skewness')

legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 6 - FRAI

figtitle = 'begin1E_PlaceFieldFRAI';
figure('Name', figtitle, 'Position', [594 467 699 420])

meanFRAI = cellfun(@mean, lapFRAI);
semFRAI = cellfun(@semfunct, lapFRAI);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:numLaps, meanFRAI(g,:), semFRAI(g,:), semFRAI(g,:), 'Color', rgb(cols{g}));
    
    leg{g} = [group(g).name ' n = ' num2str(length(lapFRAI{g})) ' fields (from the first lap)'];
end %group

xlim([0 numLaps+1])
xlabel('Laps')
xticks(1:1:numLaps)

ylabel('Firing rate aymmsetry index')
%zero_line;
ylim([-0.6 0.6])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

cd(curDir)

%end %function