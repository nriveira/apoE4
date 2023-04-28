function adprojCircTrack_x_phase_precession(group)
% function fmr1CircTrack_x_phase_precession(group)
%
% PURPOSE:
%   Plot phase precession of neurons on the circle track from WT and KO
%   groups.
%
% INPUT:
%   group struct
%
% OUTPUT:
%   Figures:
%       F(n-n) (optional): Plot showing the phase precession from the whole
%           day for each unit (only spikes included in analysis).
%       F1: Heat plots of spike counts for theta phase x normalized
%           position in place field (as in Bieri et al. 2014). Note: Spike
%           counts will be different for each group since number of cells
%           and cell firing rates are not necessarily the same.
%       F2: Slope of circular-linear regression line.
%       F3: Phase offset of circular-linear regression line.
%       F4: r^2 value of circular-linear regression line.
%
% OPTIONS:
%   See code for options.
%
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1; %1 to save figs, 0 to not

makeUnitPlots = 1; %1 to make the plots, 0 to not

downSampSpkCnts = 0; %to down sample spike counts on the first plot to be equal between groups

degBinSz = 10; %degrees per bin for theta phase
posBinSz = 0.01; %normalized distance into place field

prepForStats = 0;

%% INITIALIZE

velFilt = 1;
durCrit = 1;
runThresh = 5; %cm/s

minSpks = 50; %minimum number of spikes across the pf for each begin for cell to be included
spatCorMin = 0.7; %minimum spat coherence between first and last run sessions

numDegBins = 720/degBinSz; %two theta cycles
numPosBins = 1/posBinSz;

spatBinSz = 4;
degCmConv = (pi*100)/360; %track has 1 m diameter

degxPos = zeros(numDegBins,numPosBins,2); %initialize, by group
spkInfo = cell(2,1); %for down sampling spikes

ppSlopes = cell(2,1); %phase precession slopes by group
ppSlopesCm = cell(2,1);
ppPhaseOffset = cell(2,1);
ppPhaseRange = cell(2,1); %phase range = phase difference between the
% spike with the highest shifted phase and the spike with the lowest
% shifted phase (shifted phase is minus phase offset)
ppPhaseRangeAlt = cell(2,1); %alt phase range - pp slope * pf size
ppR2 = cell(2,1); %r2 values from circular regression

cols = {'Blue', 'Red'};
groupNames = {'apoE4', 'WT'};

saveDir = 'E:\resultsFeb2023_AD_WT\phasePrecession';
curDir = pwd;

%% GET DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        
        for d = 1:length(group(g).rat(r).day)
%         if g == 1
%             d = 5;
%         else
%             d = 3;
%         end    
%        for d = d
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                fprintf('\t\t\tUnit %d/%d\n', u, length(group(g).rat(r).day(d).xAllBeginUnitInfo))
                tetNum = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID(1);
                clustNum = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID(2);
                pf = group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf;
                if isempty(pf)
                    continue %to next cell
                end %no pf
                
                spatCorrCheck = corrcoef(group(g).rat(r).day(d).begin(1).unit(u).smRateMap, group(g).rat(r).day(d).begin(4).unit(u).smRateMap);
                if spatCorrCheck < spatCorMin 
                    continue
                end %spat corr check
                
                for p = 1:length(pf) %place field ind
                    spkCheck = 1; %initialize as a pass
                    
                    if spkCheck == 0
                        continue %to next pf
                    end %spk check
                    
                    if pf(p).radPos(end) > pf(p).radPos(1)
                        crossZero = 0;
                    else
                        crossZero = 1;
                    end %determine if pf crosses 0
                    
                    pfLen = rad2deg(circ_dist(deg2rad(pf(p).radPos(end)), deg2rad(pf(p).radPos(1))));
                    pfLenCm = pfLen * degCmConv;
                    
                    uPhisAll = [];
                    uNormPosAll = [];
                    tetNum = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID(1);
                    
                    %for b = 1:4
                    for b = 4    
                        begPasses = []; %initialize for this begin
                        begSpks = [];
                        
                        radPos = group(g).rat(r).day(d).begin(b).radPos;
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        runSpeed = smooth_runspeed(get_runspeed(coords));
                        spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                        
                        cd(group(g).rat(r).day(d).begin(b).dir)
                        lfpStruct = read_in_lfp(['CSC' num2str(tetNum) '.ncs']);
                        load(['CSC' num2str(tetNum) '_narrowThetaLfp.mat'], 'filtLfp')
                        lfpStruct.narrowThetaLfp = filtLfp;
                        
                        pfPassBnry = zeros(1,size(radPos,1));
                        if crossZero == 0
                            pfPassBnry(radPos(:,2)>= pf(p).radPos(1) & radPos(:,2)<= pf(p).radPos(end)) = 1;
                        else
                            pfPassBnry(radPos(:,2)>= pf(p).radPos(1) & radPos(:,2)<=360) = 1; %from start-360 and 0-end
                            pfPassBnry(radPos(:,2)>=0 & radPos(:,2)<= pf(p).radPos(end)) = 1;
                        end %cross zero
                        
                        pfPassChunks = bwconncomp(pfPassBnry);
                        for c = 1:length(pfPassChunks.PixelIdxList)
                            tmpInds = pfPassChunks.PixelIdxList{c};
                            
                            passDist = abs(rad2deg(circ_dist(deg2rad(radPos(tmpInds(1),2)), deg2rad(radPos(tmpInds(end),2)))));
                            
                            % Make sure rat traverses almost the whole field
                            if passDist >= pfLen - 5 % a little room for bin rounding error
                                begPasses = [begPasses; radPos(tmpInds(1),1) radPos(tmpInds(end),1)];
                            end %through whole field
                        end %chunks
                        
                        for ps = 1:size(begPasses,1)
                            
                            psSpkTms = spkTms(spkTms >= begPasses(ps,1) & spkTms <= begPasses(ps,2));
                            if isempty(psSpkTms)
                                continue
                            end
                            psRadPos = radPos(radPos(:,1)>= begPasses(ps,1) & radPos(:,1)<= begPasses(ps,2),:);
                            psCoords = coords(coords(:,1)>= begPasses(ps,1) & coords(:,1) <= begPasses(ps,2),:);
                            psSpeed = runSpeed(coords(:,1)>= begPasses(ps,1) & coords(:,1) <= begPasses(ps,2),:);
                            
                            if min(psSpeed(:,2)) < runThresh
                                continue
                            end %check run speed through the pass
                            
                            [~,~,~,spkCnts] = get_ratemap_circtrack(psSpkTms, psCoords, psRadPos, spatBinSz, velFilt, durCrit);
                            if sum(spkCnts(pf(p).inds)) < 3 %at least 3 spikes
                                continue
                            end %at least 3 spikes
                            if length(find(spkCnts(pf(p).inds))) < 2 %in at least two different position bins
                                continue
                            end %at least 2 pos bins
                            
                            begSpks = [begSpks psSpkTms'];
                        end %pass
                        if isempty(begSpks)
                            continue
                        end %no spikes
                        
                        pp = phase_precession_circtrack_pf(begSpks, radPos, coords, pf(p), lfpStruct);
                        
                        if ~isempty(pp.spkTms)
                            uPhisAll = [uPhisAll pp.spkPhis];
                            uNormPosAll = [uNormPosAll pp.normSpkPos];
                        end %if
                    end %begin
                    
                    if length(uNormPosAll) < minSpks
                        continue %to next pf
                    end %not enough spikes
                    
                    tmpDegxPos = zeros(numDegBins/2,numPosBins); %initialize %first 360 degrees
                    degBins = 0:degBinSz:360;
                    
                    allPhis = pp.spkPhis; %get all of the spike theta phases
                    
                    for db = 1:numDegBins/2
                        degRange = [degBins(db) degBins(db+1)];
                        pullPos = uNormPosAll(uPhisAll > degRange(1) & uPhisAll <= degRange(2));
                        
                        tmpPosBins = histcounts(pullPos, 0:posBinSz:1);
                        tmpDegxPos(db,:) = tmpPosBins;
                    end %degrees bins

                    tmpTwoCycle = [tmpDegxPos; tmpDegxPos];
                    degxPos(:,:,g) = degxPos(:,:,g) + tmpTwoCycle;
                    
                    beta = CircularRegression(uNormPosAll, deg2rad(uPhisAll));
                    slope = rad2deg(beta(1));
                    phaseOff = rad2deg(wrapTo2Pi(beta(2)));
                    cmSlope = slope / pfLenCm;
                    
                    altPhaseRange = cmSlope * pfLenCm;
                    
                    phaseRange = max(uPhisAll) - min(uPhisAll);
                    % R2 = circ_corrcl(deg2rad(uPhisAll), uNormPosAll);
                    tmpCal = beta(1)*uNormPosAll + beta(2);
                    R2 = circLinRegress_r2(deg2rad(uPhisAll), tmpCal);
                    if isnan(R2)
                        R2 = [];
                    end
                    
                    if makeUnitPlots == 1
                        cd(saveDir)
                        cd('unitPlots')
                        cd(group(g).name)
                        try cd(group(g).rat(r).name)
                        catch
                            mkdir(group(g).rat(r).name)
                            cd(group(g).rat(r).name)
                        end
                        figtitle = [group(g).rat(r).day(d).name '_TT' num2str(tetNum) '_' num2str(clustNum)];
                        figure('Name', figtitle)
                        
                        hold on
                        plot(uNormPosAll, uPhisAll, 'k.')
                        plot(uNormPosAll, (uPhisAll+360), 'k.')
                        ylim([0 720])
                        
                        calphase = slope*[0 1] + phaseOff;
                        plot([0 1], calphase, 'r')
                        plot([0 1], calphase+360, 'r')
                        if max(calphase+360) < 720
                            plot([0 1], calphase+720, 'r')
                        end %in range
                        
                        ylim([0 720])
                        ylabel('Theta phase (deg)')
                        
                        xlim([0 1])
                        xlabel({'Normalized position in place field', ['slope = ' num2str(round(slope)) ' deg/pf | ' num2str(round(cmSlope,1)) ' deg/cm']})
                        
                        title(['TT' num2str(tetNum) '\_' num2str(clustNum) ])
                        
                        if saveOrNot == 1
                            saveas(gcf, figtitle, 'epsc');
                            saveas(gcf, figtitle, 'fig');
                            saveas(gcf, figtitle, 'png');
                        end %save
                    end %make unit plots
                    
                    ppSlopes{g} = [ppSlopes{g} slope];
                    ppSlopesCm{g} = [ppSlopesCm{g} cmSlope];
                    ppPhaseOffset{g} = [ppPhaseOffset{g} phaseOff];
                    ppPhaseRange{g} = [ppPhaseRange{g} phaseRange];
                    ppPhaseRangeAlt{g} = [ppPhaseRangeAlt{g} altPhaseRange];
                    ppR2{g} = [ppR2{g} R2];
                end %place field ind
            end %unit
            close all
        end %day
    end %rat
end %group

cd(saveDir)
keyboard

%% FIG 1 - PHASE X POS HEATMAP

figtitle = 'AllBegin4_PhasePrecession_ThetaPhasexPostion';
figure('Name', figtitle, 'Position', [529 461 940 420])

% maxSpkCnt = max(degxPos(:));
% normDegxPos = degxPos ./ maxSpkCnt;

for g = 1:2
    normSpkCnts = degxPos(:,:,g) ./ max(max(degxPos(:,:,g)));
    
    subplot(1,2,g)
    imagesc(linspace(0, 1, size(degxPos,2)), linspace(0, 720, size(degxPos,1)), normSpkCnts)
    colormap(jet)
    axis xy
    
    xlabel('Position in place field')
    
    if g == 1
        ylabel('Theta phase (deg)')
    end
    
    cbr = colorbar;
    ylabel(cbr, 'Normalized spike Counts')
%     ax = gca;
    caxis([0 1])
    title([groupNames{g} ' (n = ' num2str(length(ppSlopes{g})) ' cells)'])
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 2 - PHASE X POS HEATMAP RECONSTRUCTED

figtitle = 'AllBegin4_PhasePrecession_ThetaPhasexPostion_reconstructed';
figure('Name', figtitle, 'Position', [558 332 748 420])

% maxSpkCnt = max(degxPos(:));
% normDegxPos = degxPos ./ maxSpkCnt;

radBinCtrs = deg2rad(degBinSz/2:degBinSz:360);
for g = 1:2
    subplot(1,2,g)
    
    com = zeros(1,size(degxPos,2)/2);
    for p = 1:length(com)
        tmpCom = circ_mean(radBinCtrs', degxPos(1:size(degxPos,1)/2,p,g));
        com(p) = wrapTo360(rad2deg(tmpCom));
    end %pos bins

    plot(linspace(0, 1, length(com)), com, 'k.')
    hold on;
    plot(linspace(0, 1, length(com)), com+360, 'k.')
    
    beta = CircularRegression(linspace(0, 1, length(com)), deg2rad(com));
    calphase = wrapTo360(rad2deg(beta(1)*[0 1] + beta(2)));
    
    plot([0 1], calphase, 'b')
    plot([0 1], calphase+360, 'b')
    
    xlabel('Position in place field')
    
    if g == 1
        ylabel('Theta phase (deg)')
    end
    
    ylim([0 720])
    
    title({[groupNames{g} ' (n = ' num2str(length(ppSlopes{g})) ' cells)'], ['slope = ' num2str(round(rad2deg(beta(1)))) ' deg/pf, phase offset = ' num2str(round(rad2deg(beta(2)))) ' deg']})
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 2 - SLOPE
% 
% figtitle = 'PhasePrecession_Slope';
% figure('Name', figtitle, 'Position', [425 436 781 420])
% 
% yLab = 'Slope (deg/pf)';
% xLabs = {group(1).name, group(2).name};
% bar_and_dotplot(ppSlopes, yLab, xLabs, cols)
% 
% subplot(1,2,2)
% ylim([-200 1000])
% if saveOrNot == 1
%     saveas(gcf, figtitle, 'epsc');
%     saveas(gcf, figtitle, 'fig');
%     saveas(gcf, figtitle, 'png');
% end

%% FIG 3 - CM SLOPE

figtitle = 'AllBegin4_PhasePrecession_Slope_cm';
figure('Name', figtitle, 'Position', [425 436 781 420])

yLab = 'Slope (deg/cm)';
xLabs = groupNames;
bar_and_dotplot(ppSlopesCm, yLab, xLabs, cols)

subplot(1,2,2)
ylim([-25 15])

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 4 - PHASE OFFSET

figtitle = 'PhasePrecession_PhaseOffset';
figure('Name', figtitle, 'Position', [425 436 781 420])

yLab = 'Phase offset (deg)';
xLabs = groupNames;

%can't use normal mean - must be circular mean

jitter = 0.1;  %for dotplot

colMatrix = []; %initialize
for c = 1:length(cols)
    colMatrix = [colMatrix; rgb(cols{c})]; %fill in for dotplot
end %cols

avgData = zeros(1,1);
semData = zeros(1,1);
for g = 1:2
    pullRadData = deg2rad(ppPhaseOffset{g});
    avgData(g) = rad2deg(circ_mean(pullRadData'));
    
    tmpStd = rad2deg(circ_std(pullRadData'));
    semData(g) = tmpStd ./ (sqrt(length(pullRadData)));
end %group

subplot(1,2,1)

bgraph = bar(avgData, 'FaceColor', 'Flat');
hold on;
errorbar(1:length(ppPhaseOffset), avgData, semData, 'Color', 'Black', 'LineStyle', 'None')
for g = 1:length(ppPhaseOffset)
    bgraph.CData(g,:) = rgb(cols{g});
end %group

ylabel(yLab)
xticklabels(xLabs)

subplot(1,2,2)
dotplot(1:length(ppPhaseOffset), ppPhaseOffset, jitter, colMatrix, repmat(rgb('Black'), length(ppPhaseOffset), 1));

for g = 1:length(ppPhaseOffset)
    line([g-.3 g+.3], [avgData(g) avgData(g)], 'Color', 'Black', 'LineWidth', 3)
end %group

ylabel(yLab)
ylim([0 360])
xticklabels(xLabs)

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 5 - PHASE RANGE - method 1

figtitle = 'PhasePrecession_PhaseRange_m1';
figure('Name', figtitle, 'Position', [425 436 781 420])

yLab = 'Phase range (deg)';
xLabs = groupNames;
bar_and_dotplot(ppPhaseRange, yLab, xLabs, cols)

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 6 - PHASE RANGE - method 2

figtitle = 'PhasePrecession_PhaseRange_m2';
figure('Name', figtitle, 'Position', [425 436 781 420])

yLab = 'Phase range (deg)';
xLabs = groupNames;
bar_and_dotplot(ppPhaseRangeAlt, yLab, xLabs, cols)

subplot(1,2,2)
ylim([-1000 1000])

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end


%% FIG 5 - R2

figtitle = 'PhasePrecession_R2';
figure('Name', figtitle, 'Position', [425 436 781 420])

yLab = 'R^2';
xLabs = groupNames;
bar_and_dotplot(ppR2, yLab, xLabs, cols)

subplot(1,2,1)
ylim([0 0.15])

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% STATS

if prepForStats == 1
    keyboard
    %Watson-Williams test, as in Bieri et al 2014, done in Matlab:
    
    wwInput = zeros(2,numPosBins);
    for g = 1:2
        pullDegxPos = degxPos(1:numDegBins/2,:,g);
        [~, maxInds] = max(pullDegxPos);
        
        degBinCtrs = degBinSz/2:degBinSz:360;
        radBinCtrs = deg2rad(degBinCtrs);
        
        wwInput(g,:) = radBinCtrs(maxInds);
        
    end %group
    
    circ_wwtest(wwInput(1,:), wwInput(2,:));
    
    statSlopes = zeros(length(ppSlopes{1})+length(ppSlopes{2}),2);
    statPhaseOffset = zeros(length(ppSlopes{1})+length(ppSlopes{2}),2);
    statPhaseRange = zeros(length(ppSlopes{1})+length(ppSlopes{2}),2);
    statR2 = zeros(length(ppSlopes{1})+length(ppSlopes{2}),2);
    
    i = 0;
    for g = 1:2
        for u = 1:length(ppSlopes{g})
            i = i + 1;
            statSlopes(i,:) = [g ppSlopes{g}(u)];
            statPhaseOffset(i,:) = [g ppPhaseOffset{g}(u)];
            statPhaseRange(i,:) = [g ppPhaseRnage{g}(u)];
            statR2(i,:) = [g ppR2{g}(u)];
        end %units
    end %group
    
    
end %prep for stats
cd(curDir)
keyboard
end %function