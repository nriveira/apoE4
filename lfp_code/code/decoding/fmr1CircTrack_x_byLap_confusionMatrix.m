function fmr1CircTrack_x_byLap_confusionMatrix(group)
% function fmr1CircTrack_x_byLap_confusionMatrix(group)
%
% PURPOSE:
%   To plot a confusion matrix a la Zheng et al., 2020, bioRxiv, Sup Fig 1b
%       This function uses the leave-one-lap-out method.
%
% INPUT:
%   group = data struct.
%
% OUTPUT:
%   F1: Day decoding error - cumulative proportion of decoding error, with
%       line representing a separate day.
%   F2: Average confusion matrices for both genotypes.
%   F3-n: Within-day confusion matrix.
%   Fn+1: Decoding error - cumulative proportion of decoding error,
%       collapsed across days, for each genotype.
%
% OPTIONS:
%   See function for options.
%
% MMD - og. JBT
% Edited/commented 1/2022
% Colgin Lab

%% OPTIONS
minFr = 1; %minimum firing rate to include cell in decoder

minCell = 20; %minimum number of simultaneously recorded cells to do the decoding

incCI = 0; %whether or not to include confidence intervals on cumulutaive error plot
nboot = 5000; % Be warned that this takes a LONG TIME. Could lower number of boot strapped samples to make it quicker

downSamp = 0; %1 to downsamp cells, 0 to not
if downSamp == 1
    newCellNum = 35; %number of cells to down sample to
end

% Plot confusion matrices for each day
%  + save options
plotConfMatsByDay = 1; %Set to 1 to plot confusion matrices for each day, with cell counts in the title
saveDayPlots = 1; %Set to 1 to save these plots automatically in...
dayPlotSaveDir = 'E:\Rat381\results\decoding_general\confusionMatricesByDay\leaveLapOutMethod';

%Options for saving the across session confusion matrix plots
saveAllSesConfMats = 1; %Set to 1 to save confusion matrix where data is concatenated across all days
allSesConfMatSaveDir = 'E:\Rat381\results\decoding_general\confusionMatrices_acrossAllSes';

% Options for saving the cumulative distribution plot
saveCumErrPlots = 1;
cumErrSaveDir = 'E:\Rat381\results\decoding_general\decodeCumulativeErrAcrossSessions';

%% INITIALIZE

% Decoding parameters
bayesWin = 40/1000; %as in Zheng et al: 40 ms sliding time window that shifted 10 ms at each step
bayesStep = 10/1000;
sampRate = 2000; %for spike raster
runThresh = 5; %cm/s

spatBinSz = 4;
velFilt = 1;
durCrit = 1;

degBinCtrs = group(1).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(degBinCtrs);
newRewLoc = 0; %used to get all plots to line up as though rewards are at 0 & 180
[~,newRewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(newRewLoc))-0));

dayDecErrFig = figure('name', 'Day Decoding Error', 'Position', [25 207.4 1220.8 506.4]);
dayDecErrorFigAll = figure('Name', 'DayDecodingError_allDays');
confMatFig = figure('name', 'Confusion Matrix', 'Position', [251         244        1421         511]);

%cumErr = cell(2,1);
cumErr = cell(1,1);

%cols = {'Blue', 'Red'};
cols = {'Blue'};
%% GET DATA

for g = 1
    fprintf('Group %d\n', g);
    
    figure(dayDecErrFig);
    subplot(1,2,g);
    xlim([0 180]);
    xlabel('Decoding Error (degrees)');
    if g == 1
        ylabel('Cumulative Proportion');
    end
    title(group(g).name);
    fix_font;
    axis square;
    hold on;
    
    allConfMats = [];
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            lapCntr = 0;
            rewLoc = group(g).rat(r).day(d).rewLocs(1);
            [~,rewInd] = min(abs(circ_dist(deg2rad(degBinCtrs), deg2rad(rewLoc))-0));
            shiftVal = newRewInd - rewInd;
            dayDecodeErr = [];
            
            % confMat = zeros(length(radBinCtrs),length(radBinCtrs)); %Confusion Matrix sums
            dayConfMat = []; %for eventually averaging across all laps in a day
            if length(group(g).rat(r).day(d).xAllBeginUnitInfo) >= minCell
                badU = [];
                uIDs = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo),2);
                
                allBegSpkCnts = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo), length(group(g).rat(r).day(d).xAllBeginUnitInfo(1).smRateMap));
                uPkPos = zeros(1,length(group(g).rat(r).day(d).xAllBeginUnitInfo));
                for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                    if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap)>= minFr %unit is bad if max firing rate in bin does not exceed 1
                        for b = 1:4
                            allBegSpkCnts(u,:) = allBegSpkCnts(u,:) + group(g).rat(r).day(d).begin(b).unit(u).spkCnts;
                        end %begin
                        uIDs(u,:) = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                        [~, pkInd] = max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap);
                        uPkPos(u) = degBinCtrs(pkInd);
                    else
                        badU = [badU u]; %#ok
                    end %bad fr
                end %units
                uIDs(badU,:) = [];
                allBegSpkCnts(badU,:) = [];
                allBegTpb = group(g).rat(r).day(d).xAllBeginTpb;
                
                for b = 1:4
                    fprintf('\t\t\t\tBegin %d/4\n', b)
                    radPos = group(g).rat(r).day(d).begin(b).radPos;
                    coords = group(g).rat(r).day(d).begin(b).coords;
                    lapTms = group(g).rat(r).day(d).begin(b).lapTms;
                    instRs = get_runspeed(coords);
                    smRs = smooth_runspeed(instRs);
                    fprintf('\t\t\t\t\tLap #...')
                    for lp = 1:size(lapTms,1)
                        lapCntr = lapCntr + 1;
                        fprintf('  %d  ', lp)
                        
                        lpCoords = coords(coords(:,1) >= lapTms(lp,1) & coords(:,1) <= lapTms(lp,2),:);
                        lpRadPos = radPos(radPos(:,1) >= lapTms(lp,1) & radPos(:,1) <= lapTms(lp,2),:);
                        
                        [~, ~, lapTpb, ~] = get_ratemap_circtrack([], lpCoords, lpRadPos, spatBinSz, velFilt, durCrit);
                        newTpb = allBegTpb - lapTpb;
                        
                        rateMaps = zeros(length(uIDs), length(group(g).rat(r).day(d).xAllBeginUnitInfo(1).smRateMap));
                        nTimeBins = round((lapTms(lp,2)-lapTms(lp,1)) * sampRate); %number of bins in spike raster
                        spkRstr = zeros(size(uIDs,1), nTimeBins);
                        uCntr = 0;
                        
                        for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                            uID = group(g).rat(r).day(d).x3BeginUnitInfo(u).ID;
                            if ismember(uID, uIDs, 'row')
                                uCntr = uCntr + 1;
                                allSpkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                                lpSpkTms = allSpkTms(allSpkTms >= lapTms(lp,1) & allSpkTms <= lapTms(lp,2),:);
                                
                                [~, ~, ~, lapSpkCnts] = get_ratemap_circtrack(lpSpkTms, lpCoords, lpRadPos, spatBinSz, velFilt, durCrit);
                                
                                newSpkCnt = allBegSpkCnts(uCntr,:) - lapSpkCnts;
                                rateMaps(uCntr,:) = newSpkCnt ./ newTpb;
                                
                                timePassed = lpSpkTms - lapTms(lp,1);
                                spkInds = round(timePassed * sampRate);
                                spkInds(spkInds<0)= []; %if spike occured before we got a time stamp for the tracking
                                spkInds(spkInds==0) = 1;
                                spkRstr(uCntr, spkInds) = 1;
                            end %unit is included
                        end %unit for new ratemaps and spike raster
                        rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
                        rateMaps = circshift(rateMaps, [0 shiftVal]);
                        lapPxn = BayesianDecoder(spkRstr,rateMaps,bayesWin,bayesStep,sampRate); %Ernie's decoder
                        
                        [nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*sampRate, bayesStep*sampRate);
                        if winStartInds(end)+bayesWin*sampRate < size(spkRstr,2)
                            nWin = nWin + 1;
                            winStartInds(end+1) = winStartInds(end)+bayesStep*sampRate; %#ok
                        end
                        winStartTms = lapTms(lp,1) + winStartInds/sampRate;
                        winEndTms = winStartTms + bayesWin;
                        % Get the rat's actual position
                        actPosn = nan(1,length(winStartTms));
                        for i = 1:length(winStartTms)
                            winSpd = mean(smRs(smRs(:,1)>=winStartTms(i) & smRs(:,1)<winEndTms(i),2));
                            if winSpd > runThresh
                                actPosn(i) = wrapTo360(rad2deg(circ_mean(deg2rad(radPos(smRs(:,1)>=winStartTms(i) & smRs(:,1)<winEndTms(i),2))))); %get mean position across this window
                            end %speed above threshold
                        end %each window
                        
                        actPosn = wrapTo360(rad2deg(circ_dist(deg2rad(actPosn), deg2rad(rewLoc))));
                        
                        lapPxn(:,nWin+1:size(lapPxn,2)) = [];
                        
                        decodedPosBins = nan(1,size(lapPxn,2));
                        decodedPosns = nan(1,size(lapPxn,2));
                        for i = 1:size(lapPxn,2)
                            if ~isnan(lapPxn(1,i))
                                tmpPpm = lapPxn(:,i);
                                [maxVal,decodedPosBins(i)] = max(tmpPpm); %max inds
                                if length(find(tmpPpm == maxVal)) > 1
                                    keyboard
                                end
                            end %not nan (spikes occured in this time bin)
                        end %time windows
                        decodedPosns(~isnan(decodedPosBins)) = degBinCtrs(decodedPosBins(~isnan(decodedPosBins)));
                        
                        decodingErr = abs(rad2deg(circ_dist(deg2rad(actPosn), deg2rad(decodedPosns))));
                        cumErr{g} = [cumErr{g} decodingErr(~isnan(decodingErr))];
                        dayDecodeErr = [dayDecodeErr decodingErr(~isnan(decodingErr))]; %#ok
                        
                        posInds = match(actPosn,degBinCtrs);
                        posInds(isnan(actPosn)) = NaN;
                        tmpConfMat = zeros(length(radBinCtrs),length(radBinCtrs));
                        posCount = zeros(length(radBinCtrs),1); %For the denominator to calc confMatrix probability
                        for p = 1:length(radBinCtrs)
                            % Add up the probability the rat is at each location for the time windows in which the rat’s true position is being indexed
                            tmpConfMat(:,p) = nansum(lapPxn(:,posInds==p),2); % Can't be zero due to precision
                            posCount(p) = posCount(p) + nansum(posInds==p);
                        end %pos bins
                        
                        tmpConfMat = tmpConfMat./repmat(posCount', length(radBinCtrs),1);
                        tmpConfMat(isinf(tmpConfMat)) = 0;
                        dayConfMat = cat(3, dayConfMat, tmpConfMat);
                        
                    end %lap
                    fprintf('\n')
                end %begin
                
                allConfMats = cat(3, allConfMats, mean(dayConfMat,3));
                try
                    [cumProp, xScale] = calc_cum_prop_ip_range(dayDecodeErr, [0 180]);
                catch; keyboard; end
                figure(dayDecErrFig);
                plot(xScale, cumProp, 'Color', [.3 .3 .3]);
                figure(dayDecErrorFigAll)
                hold on;
                plot(xScale, cumProp, 'Color', rgb(cols{g}));
                % Plot and save for individual days
                if plotConfMatsByDay ==1
                    figName = [group(g).name '_' group(g).rat(r).name '_day' group(g).rat(r).day(d).name];
                    
                    figure('name', figName);
                    colMap = define_cust_color_map('white', 'black', 200);
                    imagesc(degBinCtrs, degBinCtrs, mean(dayConfMat,3))
                    axis xy
                    axis square
                    colormap(colMap);
                    cbr = colorbar;
                    ylabel(cbr, 'Probability')
                    
                    title({[group(g).rat(r).name]; [group(g).rat(r).day(d).name]; [num2str(length(uIDs)) ' units']})
                    xlabel('Actual Position (degrees)');
                    ylabel('Decoded Position (degrees)');
                    
                    if saveDayPlots == 1
                        curDir = pwd;
                        cd(dayPlotSaveDir);
                        saveas(gcf, figName, 'epsc')
                        saveas(gcf, figName, 'png')
                        saveas(gcf, figName, 'fig')
                        cd(curDir);
                    end
                    
                end %plot each day
            end %enough cells
        end %day
    end %rat
    
    figure(confMatFig);
    %subplot(1,2,g);
    colMap = define_cust_color_map('white', 'black', 200);
    imagesc(degBinCtrs, degBinCtrs, nanmean(allConfMats,3))
    colormap(colMap);
    axis xy;
    xlabel('Actual Position (degrees)');
    ylabel('Decoded Position (degrees)');
    title(group(g).name);
    axis square;
end %group
keyboard
xlim([0 180]);
xlabel('Decoding Error (degrees)');
ylabel('Cumulative Proportion');
if saveCumErrPlots == 1
    curDir = pwd;
    cd(cumErrSaveDir)
    saveas(gcf, 'DayDecodingError_allDays', 'epsc');
    saveas(gcf, 'DayDecodingError_allDays', 'fig');
    saveas(gcf, 'DayDecodingError_allDays', 'png');
    cd(curDir)
end


figure(dayDecErrorFigAll)


figure(confMatFig);
same_caxis
cbr = colorbar;
cbr.Position = [0.9 0.1174 0.0188 0.8004];
ylabel(cbr, 'Probability')

if saveAllSesConfMats == 1
    curDir = pwd;
    cd(allSesConfMatSaveDir);
    figure(confMatFig);
    figtitle = 'confusionMatrixAcrossSessions_lapMethod';
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir);
end

if saveCumErrPlots == 1
    curDir = pwd;
    cd(cumErrSaveDir);
    figure(dayDecErrFig);
    figtitle = 'dayDecodingError_lapMethod';
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir);
end

figtitle = 'DecodingErrorAcrossSession_lapMethod';
figure;
hold on;
nboot=10;
for g = 1
    
    tmpData = cumErr{g};
    ipRange = [0 180];
    [cumProp, xScale] = calc_cum_prop_ip_range(tmpData, ipRange);
    if incCI == 1
        func = @(x) calc_cum_prop_ip_range(x,ipRange);
        ci = bootci(nboot,{func,tmpData}, 'type', 'per');
        
        line_handle(g) =  plot_filled_ci(xScale, cumProp, ci, rgb(cols{g}));
    else
        line_handle(g) = plot(xScale, cumProp, 'LineWidth', 1.5, 'Color', rgb(cols{g})); %#ok
    end %whether to include CI
    
end
legend(line_handle, {'WT', 'KO'}, 'Location', 'SouthEast');
axis square;
xlabel('Error (degrees)');
xlim([0 180])
ylabel('Cumulative Proportion');
title('Decoding Error Across the Entire Session');

if saveCumErrPlots == 1
    curDir = pwd;
    cd(cumErrSaveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end



end %function