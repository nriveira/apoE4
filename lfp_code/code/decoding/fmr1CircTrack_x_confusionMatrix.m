function fmr1CircTrack_x_confusionMatrix(group)
% function fmr1CircTrack_x_confusionMatrix(group)
%
% PURPOSE:
%   To plot a confusion matrix a la Zheng et al., 2020, bioRxiv, Sup Fig 1b
%       This function decodes across begin 4 using ratemaps from begins
%       1-3.
%
% INPUT:
%    group = project uber data struct
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

minCell = 20; %minimum number of simultaneously recorded cells to do the decoding

downSamp = 0; %1 to downsamp, 0 to not
if downSamp == 1
    newCellNum = 40; %number of cells to down sample to
end

% Plot confusion matrices for each day
%  + save options
plotConfMatsByDay = 1; %Set to 1 to plot confusion matrices for each day, with cell counts in the title
saveDayPlots = 1; %Set to 1 to save these plots automatically in...
dayPlotSaveDir = 'E:\FMR1_CIRCTRACK\RESULTS\DECODING_GENERAL\confusionMatricesByDay\leaveBeginOutMethod';

%Options for saving the across session confusion matrix plots
saveAllSesConfMats = 1; %Set to 1 to save confusion matrix where data is concatenated across all days
allSesConfMatSaveDir = 'E:\FMR1_CIRCTRACK\RESULTS\DECODING_GENERAL\confusionMatrices_acrossAllSes';

% Options for saving the cumulative distribution plot
saveCumErrPlots = 1;
cumErrSaveDir = 'E:\FMR1_CIRCTRACK\RESULTS\DECODING_GENERAL\decodeCumulativeErrAcrossSessions';

%% INITIALIZE

% Decoding parameters
bayesWin = 40/1000; %as in Zheng et al: 40 ms sliding time window that shifted 10 ms at each step
bayesStep = 10/1000;
sampRate = 2000; %for spike raster
runThresh = 5; %cm/s

radBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
newRewLoc = 0; %used to get all plots to line up as though rewards are at 0 & 180
[~,newRewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(newRewLoc))-0));

dayDecErrFig = figure('name', 'Day Decoding Error', 'Position', [25 207.4 1220.8 506.4]);
confMatFig = figure('name', 'Confusion Matrix', 'Position', [251         244        1421         511]);
tmpPlotData = [];

tmpMax = zeros(1,2);

%% GET DATA

for g = 1:2
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
            if length(group(g).rat(r).day(d).x3BeginUnitInfo) >= minCell
                
                confMat = zeros(length(radBinCtrs),length(radBinCtrs)); %Confusion Matrix sums
                posCount = zeros(length(radBinCtrs),1); %For the denominator to calc confMatrix probability
                dayDecodeErr = [];
                
                % Get ratemaps
                rateMaps = zeros(length(group(g).rat(r).day(d).x3BeginUnitInfo), length(group(g).rat(r).day(d).x3BeginUnitInfo(1).smRateMap));
                
                badU = [];
                uIDs = zeros(length(group(g).rat(r).day(d).x3BeginUnitInfo),2);
                for u = 1:length(group(g).rat(r).day(d).x3BeginUnitInfo)
                    if max(group(g).rat(r).day(d).x3BeginUnitInfo(u).rateMap)>=1
                        rateMaps(u,:) = group(g).rat(r).day(d).x3BeginUnitInfo(u).smRateMap; %Smoothed ratemap
                        uIDs(u,:) = group(g).rat(r).day(d).x3BeginUnitInfo(u).ID;
                    else
                        badU = [badU u]; %#ok
                    end
                end
                rateMaps(badU,:) = [];
                rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
                uIDs(badU,:) = [];
                
                if downSamp == 1 && size(rateMaps,1) > newCellNum
                    rs = RandStream('mlfg6331_64'); %for reproducibility
                    y = datasample(rs, 1:size(rateMaps,1) ,newCellNum, 'Replace', false); %randomly sample without replaceent
                    y = sort(y); %put it back in order
                    
                    uIDs = uIDs(y,:); %just keep these
                    rateMaps = rateMaps(y,:);
                end %down sample
                
                % Shift the rate maps so, across all data, reward locations are always plotted at 0 & 180
                %   So plots here will always look like rat goes from 0 to 360
                rewLoc = group(g).rat(r).day(d).rewLocs(1);
                [~,rewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(rewLoc))-0));
                shiftVal = newRewInd - rewInd;
                rateMaps = circshift(rateMaps, [0 shiftVal]);
                
                b = 4; %decode across begin 4, not including in ratemaps for decoder
                fprintf('\t\t\tBegin %d\n', b)
                
                %Get rat's speed
                radPos = group(g).rat(r).day(d).begin(b).radPos;
                instRs = get_runspeed(group(g).rat(r).day(d).begin(b).coords);
                smRs = smooth_runspeed(instRs);
                
                %Start and end time for this 'begin'
                startBeg = radPos(1,1);
                endBeg = radPos(end,1);
                
                %Set up spike raster for decoding
                spkRstr = zeros(size(uIDs,1), round((endBeg-startBeg)*sampRate));
                
                uCntr = 0; %can't just use u since we discard some units due to firing rate (or down sampling
                for u = 1:length(group(g).rat(r).day(d).begin(b).unit)
                    uID = group(g).rat(r).day(d).x3BeginUnitInfo(u).ID;
                    if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap (or down samp)
                        uCntr = uCntr + 1;
                        
                        allSpkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                        allSpkTms = allSpkTms(allSpkTms>startBeg & allSpkTms<endBeg); %remove any that might not fit within our tracking bounds
                        timePassed = allSpkTms - startBeg;
                        spkInds = round(timePassed * sampRate);
                        spkInds(spkInds==0) = []; %in case a spike occurred before the first video time-stamp
                        spkRstr(uCntr,spkInds) = 1;
                        
                        if ~isempty(find(diff(spkInds) == 0))
                            keyboard
                        end
                    end %%unit wasn't discarded
                end %unit
                
                % Bayesian Decoder
                ppm = BayesianDecoder(spkRstr, rateMaps, bayesWin, bayesStep, sampRate);
                ppm(isnan(ppm)) = 1 / size(rateMaps,2);
                
                % Figure out the number of time windows, their start and middle times
                [nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*sampRate, bayesStep*sampRate);
                if winStartInds(end)+bayesWin*sampRate < size(spkRstr,2)
                    nWin = nWin + 1;
                    winStartInds(end+1) = winStartInds(end)+bayesStep*sampRate; %#ok
                end
                winStartTms = startBeg + winStartInds/sampRate;
                winEndTms = winStartTms + bayesWin;
                % Get the rat's actual position
                actPosn = nan(1,length(winStartTms));
                for i = 1:length(winStartTms)
                    winSpd = mean(smRs(smRs(:,1)>=winStartTms(i) & smRs(:,1)<winEndTms(i),2));
                    if winSpd > runThresh
                        actPosn(i) = rad2deg(circ_mean(deg2rad(radPos(smRs(:,1)>=winStartTms(i) & smRs(:,1)<winEndTms(i),2)))); %get mean position across this window
                    end
                end
                
                % Shift everything so laps always show rat going from 0 to 360
                actPosn = rad2deg(circ_dist(deg2rad(actPosn), deg2rad(rewLoc)));
                actPosn(actPosn<0) = actPosn(actPosn<0) + 360;
                
                % Adjust the number of time windows of the PPM since it has the potential to be miscalculated by the decoder
                ppm(:,nWin+1:size(ppm,2)) = [];
                
                % Extract the most probable position from PPM
                decodedPosBins = nan(1,size(ppm,2));
                decodedPosns = nan(1,size(ppm,2));
                for i = 1:size(ppm,2)
                    if sum(sum(diff(ppm(:,i)))) ~= 0 %if probability wasn't equal for all bins (no spikes)
                        tmpPpm = ppm(:,i);
                        [maxVal,decodedPosBins(i)] = max(tmpPpm); %max inds
                        if length(find(tmpPpm == maxVal)) > 1
                            keyboard
                        end
                    end
                end
                decodedPosns(~isnan(decodedPosBins)) = radBinCtrs(decodedPosBins(~isnan(decodedPosBins)));
                
                %Calculate error of decoding
                decodingErr = abs(rad2deg(circ_dist(deg2rad(actPosn), deg2rad(decodedPosns))));
                dayDecodeErr = [dayDecodeErr decodingErr(~isnan(decodingErr))]; %#ok
                
                %Convert actual positions to the closest index
                posInds = match(actPosn,radBinCtrs);
                posInds(isnan(actPosn)) = NaN;
                
                for p = 1:length(radBinCtrs)
                    % Add up the probability the rat is at each location for the time windows in which the rat’s true position is being indexed
                    confMat(:,p) = confMat(:,p) + nansum(ppm(:,posInds==p),2); % Can't be zero due to precision
                    %The number of time windows for which the rat's true position was equal to the indexed position
                    posCount(p) = posCount(p) + nansum(posInds==p);
                end
                %**NOTE: If there isn't a single spike, the location at which the rat spent the most time will have the greatest values in the confusion matrix
                %          - The decoder will say probability is equal for all locations if no spiking occurred during a time window.
                %          - If we then add up all those time windows split by where the rat actually was, the spot that he occupied the most
                %            will have the highest values.
                
                
                % Plot cumulative proportion of error for each session
                [cumProp, xScale] = calc_cum_prop_ip_range(dayDecodeErr, [0 180]);
                figure(dayDecErrFig);
                plot(xScale, cumProp, 'Color', [.3 .3 .3]);
                if isempty(tmpPlotData)
                    tmpPlotData = cumProp;
                else
                    tmpPlotData = [tmpPlotData; cumProp]; %#ok
                end
                
                % Calculate the confusion matrix for this day
                confMat = confMat./repmat(posCount', length(radBinCtrs),1);
                confMat(isinf(confMat)) = 0;
                
                % Concatenate across days
                allConfMats = cat(3, allConfMats, confMat);
                
                % Plot and save for individual days
                if plotConfMatsByDay ==1
                    figName = [group(g).name '_' group(g).rat(r).name '_day' group(g).rat(r).day(d).name];
                    
                    if downSamp == 1
                        figName = [group(g).name '_' group(g).rat(r).name '_day' group(g).rat(r).day(d).name '_downSamp' num2str(newCellNum) 'units'];
                    end %if
                    
                    figure('name', figName);
                    colMap = define_cust_color_map('white', 'black', 200);
                    imagesc(radBinCtrs, radBinCtrs, confMat)
                    colormap(colMap);
                    %                     caxis([0 .35]);
                    %                 title({['Day ' group(g).rat(r).day(d).name]; ['#units = ' num2str(size(spkRstr,1))]})
                    fprintf('\t\t\tUnit check: %d units\n', size(spkRstr,1))
                    %                 title({[group(g).rat(r).name]; [group(g).rat(r).day(d).name]})
                    if isfield(group(g).rat(r).day(d), 'dsNUnit')
                        title({[group(g).rat(r).name]; [group(g).rat(r).day(d).name]; ['down-sampled to ' num2str(group(g).rat(r).day(d).dsNUnit) ' units']})
                    else
                        title({[group(g).rat(r).name]; [group(g).rat(r).day(d).name]; [num2str(length(uIDs)) ' units']})
                    end %if
                    axis xy;
                    xlabel('Actual Position (degrees)');
                    ylabel('Decoded Position (degrees)');
                    axis square;
                    colorbar;
                    if saveDayPlots == 1
                        curDir = pwd;
                        cd(dayPlotSaveDir);
                        %                     print(figName, '-dpng');
                        saveas(gcf, figName, 'epsc')
                        saveas(gcf, figName, 'png')
                        saveas(gcf, figName, 'fig')
                        cd(curDir);
                    end
                end
                
                
            end %if there are enough cells to decode this day
        end %day
        
    end %rat
    
    
    figure(confMatFig);
    subplot(1,2,g);
    colMap = define_cust_color_map('white', 'black', 30);
    imagesc(radBinCtrs, radBinCtrs, nanmean(allConfMats,3))
    colormap(colMap);
    axis xy;
    xlabel('Actual Position (degrees)');
    ylabel('Decoded Position (degrees)');
    title(group(g).name);
    axis square;
end %group

figure(confMatFig);
same_caxis
cbr = colorbar;
cbr.Position = [0.9 0.1174 0.0188 0.8004];
ylabel(cbr, 'Probability')

if saveAllSesConfMats == 1
    curDir = pwd;
    cd(allSesConfMatSaveDir);
    figure(confMatFig);
    figtitle = 'confusionMatrixAcrossSessions';
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir);
end

if saveCumErrPlots == 1
    curDir = pwd;
    cd(cumErrSaveDir);
    figure(dayDecErrFig);
    figtitle = 'dayDecodingError';
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir);
end

figure;
hold on;
for g = 1:2
    figLine(g) = plot(xScale, tmpPlotData(g,:), 'LineWidth', 1.5); %#ok
end
legend(figLine, {'WT', 'KO'}, 'Location', 'SouthEast');
axis square;
xlabel('Error (degrees)');
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



end %fnctn