function adprojCircTrack_x_lapByLap_phase_precession(group)
% function fmr1CircTrack_x_lapByLap_phase_precession(group)
%
% PURPOSE:
%   Make plots to compare lap-by-lap phase precession measures.
%
% INPUT:
%   group struct
%
% OUTPUT:
%   F(n-n)(optional): Plot for each unit comparing lap to overall average from the
%       first 10 (or selected number of laps).
%   F1: Slope of the phase precession using two unit measurements:
%       - deg/place field
%       - deg/cm on the track
%   F2: Phase offset across laps.
%   F3: Phase range using 2 methods 
%       - Phase difference between the spike withthe highest shifted phase
%       and the spike with the lowest shifted phase (Feng et al. 2015).
%       - Multiplication of phase precession slope and the size of the
%       single-lap place field (Schmidt et al. 2009).
%   F4: R2 value of phase precession using circular-linear regression
%       across laps.
%
% OPTIONS:
%   See function for options.
% 
% MMD
% 8/2021
% Colgin Lab

%% OPTIONS

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

%numLaps = minLapDay - 1; %number of laps to plot
numLaps = 20; %number of laps to plot
saveDir = 'E:\resultsFeb2023_AD_WT\phasePrecession\byLap'; %for saving the pahse precession figs
% pfSaveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\placeCellProperties\byLap';
saveOrNot = 1;

makeUnitPlots = 1;

prepForStats = 0;

%% INITIALIZE

spatBinSz = 4;
velFilt = 1;
durCrit = 1;

bound = 2;

lapSlopes = cell(2,numLaps); %group x first 10 laps individually
bSlopes = cell(2,1); %all together
cmLapSlopes = cell(2,numLaps);
cmbSlopes = cell(2,1);
lapPhaseOff = cell(2,numLaps);
bPhaseOff = cell(2,1);
lapPhaseRange = cell(2,numLaps);
bPhaseRange = cell(2,1);
altLapPhaseRange = cell(2,numLaps);
altBPhaseRange = cell(2,1);
lapR2 = cell(2,numLaps);
bR2 = cell(2,1);

% lapFirRate = cell(2,numLaps);
% bFirRate = cell(2,1);
% lapPkFir = cell(2,numLaps);
% bPkFir = cell(2,1);
% lapPfCent = cell(2,numLaps);
% lapPfSizeRatio = cell(2,numLaps);
% lapPfSkew = cell(2,numLaps); 
% lapFRAI = cell(2,numLaps);

degCmConv = (pi*100)/360; %track has 1 m diameter

cols = {'Blue', 'Red'};

if prepForStats == 1
    statSlopes = [];
    statSlopesCm = [];
    statPhaseOff = [];
    statPhaseRange = [];
    statAltPhaseRange = [];
    statR2 = [];
    
    statFr = [];
    statPkPf = [];
    statPfCent = [];
    statPfSizeRat = [];
    statPfSkew = [];
    statFRAI = [];
end %initialize stats

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
                
                if ~isempty(xBegPf)
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
                            
                            for ps = 1:size(pfPasses{b},1) %stable = fire 3 spike in 2 dif bins on each lap
                                passCntr = passCntr + 1;
                                psSpkTms = spkTms(spkTms>= pfPasses{b}(ps,1) & spkTms<=pfPasses{b}(ps,2));
                                psRadPos = radPos(radPos(:,1)>=pfPasses{b}(ps,1) & radPos(:,1)<=pfPasses{b}(ps,2),:);
                                psCoords = coords(coords(:,1)>=pfPasses{b}(ps,1) & coords(:,1)<=pfPasses{b}(ps,2),:);
                                [tmpMap,~,~,spkCnts] = get_ratemap_circtrack(psSpkTms, psCoords, psRadPos, spatBinSz, velFilt, durCrit);
                                
                                if sum(spkCnts(xBegPf(p).inds)) < 3 %at least 3 spikes
                                    stableCheck = 0;
                                    break %out of pass loop
                                end %at least 3 spikes
                                if length(find(spkCnts(xBegPf(p).inds))) < 2 %in at least two different position bins
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
                            end
                        end %begin
                        
                        if stableCheck == 1
                            if prepForStats == 1
                                pfStatSlope = [g nan(1,numLaps)]; %initialize
                                pfCmStatSlope = [g nan(1,numLaps)];
                                pfStatPO = [g nan(1,numLaps)];
                                pfStatPR = [g nan(1,numLaps)];
                                pfStatAltPR = [g nan(1,numLaps)];
                                pfStatR2 = [g nan(1,numLaps)];
                                
%                                 pfStatFr = [g nan(1,numLaps)];
%                                 pfStatPkFr = [g nan(1,numLaps)];
%                                 pfStatCent = [g nan(1,numLaps)];
%                                 pfStatSizeRat = [g nan(1,numLaps)];
%                                 pfStatSkew = [g nan(1,numLaps)];
%                                 pfStatFRAI = [g nan(1,numLaps)];
                            end %stats
                            
                            allSpkPhis = []; %for getting all for entire day plot later
                            allSpkPos = [];
                            allSpkPosCm = [];
                            passCntr = 0; %re-initialize
                            if makeUnitPlots == 1
                                figtitle = [group(g).rat(r).day(d).name '_TT' num2str(tetNum) '_' num2str(clustNum)];
                              figure('Name', figtitle, 'Position', [680 114 747 864]);
                            end %make fig for unit plots
                            for b = 1:4
                                radPos = group(g).rat(r).day(d).begin(b).radPos;
                                coords = group(g).rat(r).day(d).begin(b).coords;
                                spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                                
                                cd(group(g).rat(r).day(d).begin(b).dir)
                                lfpStruct = read_in_lfp(['CSC' num2str(tetNum) '.ncs']);
                                load(['CSC' num2str(tetNum) '_narrowThetaLfp.mat'])
                                lfpStruct.narrowThetaLfp = filtLfp; %get our theta
                                
                                for ps = 1:size(pfPasses{b},1)
                                    passCntr = passCntr + 1;
                                    cutSpkTms = spkTms(spkTms >= pfPasses{b}(ps,1) & spkTms <=pfPasses{b}(ps,2)); %pull out spike time just from this pass
                                    
                                    pp = phase_precession_circtrack_pf(cutSpkTms, radPos, coords, xBegPf(p), lfpStruct);
                                    
                                    allSpkPhis = [allSpkPhis pp.spkPhis];
                                    allSpkPos = [allSpkPos pp.normSpkPos];
                                    allSpkPosCm = [allSpkPosCm pp.spkPos];
                                    
                                    if passCntr <= numLaps %store
                                        lapSlopes{g,passCntr} = [lapSlopes{g,passCntr} rad2deg(pp.stats(1))];
                                        calphase = pp.stats(1) * (0:1) + pp.stats(2);
                                        cmSlope = (calphase(2) - calphase(1)) / pfLenCm;
                                        cmLapSlopes{g,passCntr} = [cmLapSlopes{g,passCntr} rad2deg(cmSlope)];
                                        
                                        lapPhaseOff{g,passCntr} = [lapPhaseOff{g,passCntr} rad2deg(pp.stats(2))];
                                        lapPhaseRange{g,passCntr} = [lapPhaseRange{g,passCntr} max(pp.spkPhis)-min(pp.spkPhis)];
                                        lapR2{g,passCntr} = [lapR2{g,passCntr} pp.stats(3)];
                                        
%                                         lapFirRate{g,passCntr} = [lapFirRate{g,passCntr} mean(rMapByPass(passCntr,:))];
%                                         lapPkFir{g,passCntr} = [lapPkFir{g,passCntr} max(rMapByPass(passCntr,:))];
%                                         
%                                         %all the pf things
                                        psSpkBins = find(rMapByPass(passCntr,:));
                                        psPfLen = abs(rad2deg(circ_dist(deg2rad(xBegPf(p).radPos(psSpkBins(1))), deg2rad(xBegPf(p).radPos(psSpkBins(end))))));
                                        psPfCent = wrapTo360(xBegPf(p).radPos(psSpkBins(1)) + psPfLen/2);
                                        tmpDist = rad2deg(circ_dist(deg2rad(pfCent), deg2rad(psPfCent)));
%                                         
%                                         pfSkew = skewness(rMapByPass(passCntr,:))/std(rMapByPass(passCntr,:))^3;
%                                         
%                                         halfSpks = ceil(length(pp.spkTms)/2);
%                                         FR1 = halfSpks/(pp.spkTms(halfSpks) - pp.spkTms(1));
%                                         FR2 = halfSpks/(pp.spkTms(end) - pp.spkTms(halfSpks));
%                                         psFRAI = (FR1-FR2)/(FR1+FR2);
%                                         
%                                         lapPfCent{g,passCntr} = [lapPfCent{g,passCntr} tmpDist*degCmConv]; %in cm
%                                         lapPfSizeRatio{g,passCntr} = [lapPfSizeRatio{g,passCntr} psPfLen/pfLen];
%                                         lapPfSkew{g,passCntr} =[lapPfSkew{g,passCntr} pfSkew];
%                                         lapFRAI{g,passCntr} = [lapFRAI{g,passCntr} psFRAI];
%                                         if isnan(psFRAI)
%                                             keyboard
%                                         end
                                        
                                        altLapPhaseRange{g,passCntr} = [altLapPhaseRange{g,passCntr} wrapTo360(rad2deg(cmSlope)*psPfLen*degCmConv)];
                                        
                                        if prepForStats == 1
                                            pfStatSlope(passCntr+1) = rad2deg(pp.stats(1)); %
                                            pfCmStatSlope(passCntr+1) = rad2deg(cmSlope);
                                            pfStatPO(passCntr+1) = rad2deg(pp.stats(2));
                                            pfStatPR(passCntr+1) = max(pp.spkPhis)-min(pp.spkPhis);
                                            pfStatR2(passCntr+1) = pp.stats(3);
                                            
%                                             pfStatFr(passCntr+1) = mean(rMapByPass(passCntr,:));
%                                             pfStatPkFr(passCntr+1) = max(rMapByPass(passCntr,:));                                            
%                                             pfStatCent(passCntr+1) = tmpDist*degCmConv;
%                                             pfStatSizeRat(passCntr+1) = psPfLen/pfLen;
%                                             pfStatSkew(passCntr+1) = pfSkew;
%                                             pfStatFRAI(passCntr+1) = psFRAI;
                                            pfStatAltPR(passCntr+1) = wrapTo360(rad2deg(cmSlope)*psPfLen*degCmConv);
                                        end %stats
                                        
                                    end %passCntr check
                                    
                                    if makeUnitPlots == 1 && passCntr == 1
                                        
                                        xVals = [min(pp.normSpkPos) max(pp.normSpkPos)];
                                        
                                        calphase = pp.stats(1) * xVals + pp.stats(2);
                                        while calphase(1) < 0
                                            calphase = calphase + 2*pi;
                                        end
                                        while calphase(1) >= 360
                                            calphase = calphase -2*pi;
                                        end %get claphase in range
                                        calphase = rad2deg(calphase);
                                        
                                        for spt = 1:2
                                            
                                            subplot(2,2,spt)
                                            if spt == 1
                                                posData = pp.normSpkPos;
                                                ttl = ['First lap - slope = ' num2str(round(rad2deg(pp.stats(1)),1)) ' deg/pf'];
                                                xLab = 'Normalized position in place field';
                                            else
                                                posData = pp.spkPos;
                                                ttl = ['First lap - slope = ' num2str(round(rad2deg(cmSlope),1)) ' deg/cm'];
                                                xLab = 'Position (cm)';
                                            end %subplot type
                                            
                                            if spt == 2 && crossZero == 1
                                                posData(posData<pfRadPos(1)*degCmConv) = posData(posData<pfRadPos(1)*degCmConv) + 360*degCmConv;
                                                plot(posData, pp.spkPhis, 'k.')
                                                hold on;
                                                xVals = [min(posData) max(posData)];
                                                line([pfRadPos(1)*degCmConv pfRadPos(1)*degCmConv], [0 720], 'LineStyle', '--', 'Color', 'k')
                                                line([(pfRadPos(end)+360)*degCmConv (pfRadPos(end)+360)*degCmConv], [0 720], 'LineStyle', '--', 'Color', 'k')
                                                ax = gca;
                                                curXticks = ax.XLim;
                                                xticks(curXticks(1):20:curXticks(2))
                                                xticklabels([curXticks(1):20:300 6:20:curXticks(2)-314])
                                            else
                                                plot(posData, pp.spkPhis, 'k.')
                                                hold on;
                                                xVals = [min(posData) max(posData)];
                                            end%not special cross zero
                                            if spt == 2 && crossZero == 0
                                                line([pfRadPos(1)*degCmConv pfRadPos(1)*degCmConv], [0 720], 'LineStyle', '--', 'Color', 'k')
                                                line([pfRadPos(end)*degCmConv pfRadPos(end)*degCmConv], [0 720], 'LineStyle', '--', 'Color', 'k')
                                                xlim([pfRadPos(1)*degCmConv - 10 pfRadPos(end)*degCmConv+10])
                                            end %add lines to spt 2
                                            if spt == 1
                                                xlim([0 1])
                                            end
                                            
                                            plot(xVals, calphase, 'r')
                                            plot(posData, pp.spkPhis +360, 'k.')
                                            plot(xVals, calphase + 360, 'r')
                                            if max(calphase+360) < 720
                                                plot(xVals, calphase+720, 'r')
                                            end %in range
                                            ylim([0 720])
                                            ylabel('Theta phase (deg)')
                                            xlabel(xLab)
                                            title(ttl)
                                        end %subplot
                                        
                                    end %make unit plots
                                end %passes
                                
                            end %begins
                            if prepForStats == 1
                                
                                statSlopes = [statSlopes; pfStatSlope];
                                statSlopesCm = [statSlopesCm pfCmStatSlope];
                                statPhaseOff = [statPhaseOff; pfStatPO];
                                statPhaseRange = [statPhaseRange; pfStatPR];
                                statR2 = [statR2; pfStatR2];
                                statAltPhaseRange = [statAltPhaseRange; pfStatAltPR];
                                
%                                 statFr = [statFr; pfStatFr];
%                                 statPkPf = [statPkPf; pfStatPkFr];
%                                 statPfCent = [statPfCent; pfStatCent];
%                                 statPfSizeRat = [statPfSizeRat; pfStatSizeRat];
%                                 statPfSkew = [statPfSkew; pfStatSkew];
%                                 statFRAI = [statFRAI; pfStatFRAI];  
                                
                            end %stats
                            
                            para = circ_lin_regress(allSpkPos, deg2rad(allSpkPhis), 2); %using CZ code and methods
                            if para(1,1) == bound || para(1,1) == - bound
                                keyboard
                            end
%                             [beta,R2,p] =
%                             CircularRegression(allSpkPos,deg2rad(allSpkPhis)); %same answer as using bound = 2
                            
                            calphase = 2*pi*para(1,1)*(0:1)+para(1,2);
                            slope = calphase(2) - calphase(1);
                            cmSlope = (calphase(2) - calphase(1)) / pfLenCm;
                            
                            if sum(calphase < 0) < 0 && sum(calphase>=2*pi) > 0
                                keyboard
                            end
                            
                            if calphase(1) < 0 %phase off between 0 and 2*pi
                                calphase = calphase + 2*pi;
                            elseif calphase(1) >= 2*pi
                                calphase = calphase -2*pi;
                            end
                            phaseOff = calphase(1);
                            R2 = circ_corrcl(deg2rad(allSpkPhis), allSpkPos);
                            
                            bSlopes{g} = [bSlopes{g} rad2deg(slope)];
                            cmbSlopes{g} = [cmbSlopes{g} rad2deg(cmSlope)];
                            bPhaseOff{g} = [bPhaseOff{g} rad2deg(phaseOff)];
                            bPhaseRange{g} = [bPhaseRange{g} max(allSpkPhis)-min(allSpkPhis)];
                            altBPhaseRange{g} = [altBPhaseRange{g} wrapTo360(rad2deg(cmSlope)*pfLen*degCmConv)];
                            bR2{g} = [bR2{g} R2];
                            
%                             bFirRate{g} = [bFirRate{g} mean(xBegRateMap(xBegPf(p).inds))];
%                             bPkFir{g} = [bPkFir{g} max(xBegRateMap(xBegPf(p).inds))];
                            
                            if makeUnitPlots == 1
                                xVals = [min(allSpkPos) max(allSpkPos)];
                                calphase = slope * xVals + phaseOff;
                                while calphase(1) < 0
                                    calphase = calphase + 2*pi;
                                end
                                while calphase(1) >= 360
                                    calphase = calphase -2*pi;
                                end %get claphase in range
                                calphase = rad2deg(calphase);
                                
                                for spt = 1:2
                                    subplot(2,2,spt+2)
                                    
                                    if spt == 1
                                        posData = allSpkPos;
                                        ttl = ['Entire day - slope = ' num2str(round(rad2deg(slope),1)) ' deg/pf'];
                                        xLab = 'Normalized position in place field';
                                    else
                                        posData = allSpkPosCm;
                                        ttl = ['Entire day - slope = ' num2str(round(rad2deg(cmSlope),1)) ' deg/cm'];
                                        xLab = 'Position (cm)';
                                    end
                                    
                                    if spt == 2 && crossZero == 1
                                        posData(posData<pfRadPos(1)*degCmConv) = posData(posData<pfRadPos(1)*degCmConv) + 360*degCmConv;
                                        plot(posData, allSpkPhis, 'k.')
                                        hold on;
                                        xVals = [min(posData) max(posData)];
                                        line([pfRadPos(1)*degCmConv pfRadPos(1)*degCmConv], [0 720], 'LineStyle', '--', 'Color', 'k')
                                        line([(pfRadPos(end)+360)*degCmConv (pfRadPos(end)+360)*degCmConv], [0 720], 'LineStyle', '--', 'Color', 'k')
                                        ax = gca;
                                        curXticks = ax.XLim;
                                        xticks(curXticks(1):20:curXticks(2))
                                        xticklabels([curXticks(1):20:300 6:20:curXticks(2)-314])
                                    else
                                        plot(posData, allSpkPhis, 'k.')
                                        hold on;
                                        xVals = [min(posData) max(posData)];
                                    end %whether cross zero case
                                    if spt == 2 && crossZero == 0
                                        line([pfRadPos(1)*degCmConv pfRadPos(1)*degCmConv], [0 720], 'LineStyle', '--', 'Color', 'k')
                                        line([pfRadPos(end)*degCmConv pfRadPos(end)*degCmConv], [0 720], 'LineStyle', '--', 'Color', 'k')
                                        xlim([pfRadPos(1)*degCmConv - 10 pfRadPos(end)*degCmConv+10])
                                    end %add lines to spt 2
                                    if spt == 1
                                        xlim([0 1])
                                    end
                                    
                                    plot(xVals, calphase, 'r')
                                    plot(posData, allSpkPhis +360, 'k.')
                                    plot(xVals, calphase + 360, 'r')
                                    if max(calphase+360) < 720
                                        plot(xVals, calphase+720, 'r')
                                    end %in range
                                    ylim([0 720])
                                    ylabel('Theta phase (deg)')
                                    xlabel(xLab)
                                    title(ttl)
                                end %subplot type
                                if saveOrNot == 1
                                    cd(saveDir)
                                    cd('unitPlots')
                                    cd(group(g).name)
                                    try
                                        cd(group(g).rat(r).name)
                                    catch
                                        mkdir(group(g).rat(r).name)
                                        cd(group(g).rat(r).name)
                                    end
                                    saveas(gcf, figtitle, 'epsc');
                                    saveas(gcf, figtitle, 'fig');
                                    try
                                        saveas(gcf, figtitle, 'png');
                                    catch
                                        keyboard
                                    end
                                end %saveOrNot
                                
                            end %add to unit plot
                        end %pf is stable throughout begin 1
                    end %pfs
                end %there is a pf for this unit
            end %unit
            close all
        end %day
    end %rat
end %group

cd(saveDir)
keyboard
%% FIG 1 - SLOPES
figtitle = 'PhasePrecessionSlope_acrossLaps';

figure('Name', figtitle, 'Position', [248 486 1408 420])

for spt = 1
    if spt == 1
        tmpData = lapSlopes;
        tmpBData = bSlopes;
        yLab = 'Slope (deg/pf)';
    else
        tmpData = cmLapSlopes;
        tmpBData = cmbSlopes;
        yLab = 'Slope (deg/cm)';
    end
    %subplot(1,2,spt)
    meanSlope = cellfun(@mean, tmpData);
    meanBeg = cellfun(@mean, tmpBData);
    semSlopes = cellfun(@semfunct, tmpData);
    semBeg = cellfun(@semfunct, tmpBData);
    
    lh = nan(1,2);
    leg = cell(1,2);
    for g = 1:2
        hold on;
        errorbar(1:numLaps, meanSlope(g,:), semSlopes(g,:), semSlopes(g,:), 'Color', rgb(cols{g}))
        
        lh(g) = errorbar(numLaps+2, meanBeg(g), semBeg(g), semBeg(g), 'Color', rgb(cols{g}));
        leg{g} = [group(g).name ' n = ' num2str(length(lapSlopes{g})) ' cells'];
    end %group
    
    xlim([0 numLaps+3])
    xlabel('Laps')
    xticks([1:numLaps numLaps+2])
    xticklabels({[1:numLaps], 'Entire day'})
    
    ylabel(yLab)
    zero_line;
    legend(lh, leg, 'Location', 'northeastoutside')
end

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 2 - PHASE OFFSET

figtitle = 'PhasePrecessionPhaseOffset_acrossLaps';

figure('Name', figtitle, 'Position', [594 467 699 420])

meanPhaseOff = cellfun(@mean, lapPhaseOff);
meanBegPO = cellfun(@mean, bPhaseOff);
semPhaseOff = cellfun(@semfunct, lapPhaseOff);
semBegPO = cellfun(@semfunct, bPhaseOff);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    errorbar(1:numLaps, meanPhaseOff(g,:), semPhaseOff(g,:), semPhaseOff(g,:), 'Color', rgb(cols{g}))
    
    lh(g) = errorbar(numLaps + 2, meanBegPO(g), semBegPO(g), semBegPO(g), 'Color', rgb(cols{g}));
    leg{g} = [group(g).name ' n = ' num2str(length(lapPhaseOff{g})) ' cells'];
end %group

xlim([0 numLaps+3])
xlabel('Laps')
xticks([1:numLaps numLaps+2])
xticklabels({[1:numLaps], 'Entire day'})

ylabel('Phase offset (deg)')
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 3 - PHASE RANGE

figtitle = 'PhasePrecessionPhaseRange_acrossLaps';

figure('Name', figtitle, 'Position', [248 486 1408 420])

ttls = {'Method 1 (Feng et al. 2015)', 'Method 2 (Schmidt et al. 2009)'};

for spt = 1:2
    if spt == 1
        lapData = lapPhaseRange;
        bData = bPhaseRange;
    else
        lapData = altLapPhaseRange;
        bData = altBPhaseRange;
    end %subplot type get data
    
    subplot(1,2,spt)
    
    meanPhaseR = cellfun(@mean, lapData);
    meanBegR = cellfun(@mean, bData);
    semPhaseR = cellfun(@semfunct, lapData);
    semBegR = cellfun(@semfunct, bData);
    
    lh = nan(1,2);
    leg = cell(1,2);
    for g = 1:2
        hold on;
        errorbar(1:numLaps, meanPhaseR(g,:), semPhaseR(g,:), semPhaseR(g,:), 'Color', rgb(cols{g}))
        
        lh(g) = errorbar(numLaps+2, meanBegR(g), semBegR(g), semBegR(g), 'Color', rgb(cols{g}));
        leg{g} = [group(g).name ' n = ' num2str(length(lapPhaseRange{g})) ' cells'];
    end %group
    
    xlim([0 numLaps+3])
    xlabel('Laps')
    xticks([1:numLaps numLaps+2])
    xticklabels({[1:numLaps], 'Entire day'})
    
    ylabel('Phase range (deg)')
%     ylim([0 360])
    legend(lh, leg, 'Location', 'northeastoutside')
    title(ttls{spt})
end %subplot type

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 4 - R2

figtitle = 'PhasePrecessionR2_acrossLaps';

figure('Name', figtitle, 'Position', [594 467 699 420])

meanR2= cellfun(@nanmean, lapR2);
meanBegR2 = cellfun(@mean, bR2);
semR2 = cellfun(@nansemfunct, lapR2);
semBegR2 = cellfun(@semfunct, bR2);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    errorbar(1:numLaps, meanR2(g,:), semR2(g,:), semR2(g,:), 'Color', rgb(cols{g}))
    
    lh(g) = errorbar(numLaps+2, meanBegR2(g), semBegR2(g), semBegR2(g), 'Color', rgb(cols{g}));
    leg{g} = [group(g).name ' n = ' num2str(length(lapR2{g})) ' cells'];
end %group

xlim([0 numLaps+3])
xlabel('Laps')
xticks([1:numLaps numLaps+2])
xticklabels({[1:numLaps], 'Entire day'})

ylabel('R^2')
ylim([0 0.6])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% STATS

if prepForStats == 1
    keyboard
end %stats

cd(curDir)

end %function