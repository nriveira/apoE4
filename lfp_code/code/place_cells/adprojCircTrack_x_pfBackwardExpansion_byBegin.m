function adprojCircTrack_x_pfBackwardExpansion_byBegin(group)
% function fmr1CircTrack_x_pfBackwardExpansion_byBegin(group)
%
% PURPOSE:
%   The purpose of this function is to look at the backwards expansion of
%   place cells in WT and KO rats. This function looks at the place field
%   across the four begins in a day.
% 
% INPUT:
%   group = data struct
% 
% OUTPUT:
%   F1: In-field firing rate and peak firing rate in place field across
%       begins.
%   F2: Place field center difference across begins.
%   F3: Peak firing position difference across begins.
%   F4: Place field skewness across begins.
% 
% MMD
% 2/2022
% Colgin Lab

%% OPTIONS

saveOrNot = 1;

minSpkPerBeg = 50; %min pf spks across the begin for cell to be included


%% INITIALIZE

saveDir = 'E:\resultsFeb2023_AD_WT\backwardExpansion';

firRate = cell(2,4); %group x begin
pkFirRate = cell(2,4);
pkPos = cell(2,4);
pfSize = cell(2,4);
pfSizeRatio = cell(2,4);
pfCent = cell(2,4);
pfSkew = cell(2,4);
pfFRAI = cell(2,4);

cols = {'Blue', 'Red'};

curDir = pwd;
cd(saveDir)

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                xBegPf = group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf;
                if isempty(xBegPf)
                    continue %to next unit
                end %no pf
                for p = 1:length(xBegPf)
                    pfInds = xBegPf(p).inds;
                    pfPos = xBegPf(p).radPos;
                    
                    stableCheck = 1; %inialize as pass
                    for b = 1:4
                        numSpks = sum(group(g).rat(r).day(d).begin(b).unit(u).spkCnts(pfInds));
                        if numSpks < minSpkPerBeg
                            stableCheck = 0;
                            break
                        end %if not enough spikes
                    end %begin
                    
                    if stableCheck == 0
                        continue %to next pf
                    end %stable check
                    
                    pfRm = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap(pfInds);
                    pfLen = rad2deg(circ_dist(deg2rad(pfPos(end)), deg2rad(pfPos(1))));
                    pfPkPos = xBegPf(p).pkPos;
                    
                    props = regionprops(true(size(pfRm)),  pfRm, 'WeightedCentroid'); %get center of mass
                    pfCentNorm = props.WeightedCentroid(1) / length(pfRm); %normalize by the size of the pf
                    xBegCent = wrapTo360(pfPos(1) + pfLen*pfCentNorm);
                    
                    for b = 1:4
                        begPfRm = group(g).rat(r).day(d).begin(b).unit(u).smRateMap(pfInds);
                        
                        begFr = mean(begPfRm);
                        [begPkFr, begPkPosInd] = max(begPfRm);
                        begPkPos = pfPos(begPkPosInd);
                        begPkPosDiff = rad2deg(circ_dist(deg2rad(begPkPos), deg2rad(pfPkPos)));
                        
                        props = regionprops(true(size(begPfRm)),  begPfRm, 'WeightedCentroid'); %get center of mass
                        begCentNorm = props.WeightedCentroid(1) / length(begPfRm); %normalize by the size of the pf
                        begCent = wrapTo360(pfPos(1) + pfLen*begCentNorm);
                        begCentDiff = rad2deg(circ_dist(deg2rad(begCent), deg2rad(xBegCent)));
                        
                        begSkew = skewness(begPfRm);
                        
                        firRate{g,b} = [firRate{g,b} begFr];
                        pkFirRate{g,b} = [pkFirRate{g,b} begPkFr];
                        pkPos{g,b} = [pkPos{g,b} begPkPosDiff];
                        pfCent{g,b} = [pfCent{g,b} begCentDiff];
                        pfSkew{g,b} = [pfSkew{g,b} begSkew];
                    end %begin
                end %place field
            end %unit
        end %day
    end %rat
end %group

cd(saveDir)
keyboard

%% FIG 1 - FIRING RATE

figtitle = 'PlaceCellFirRate_byBeg';
figure('Name', figtitle, 'Position', [248 486 1408 420])

subplot(1,2,1)
meanData = cellfun(@mean, firRate);
semData = cellfun(@semfunct, firRate);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:4, meanData(g,:), semData(g,:), semData(g,:), 'Color', rgb(cols{g}));
    leg{g} = [group(g).name ' n = ' num2str(length(firRate{g})) ' cells'];
end %group

xlim([0 5])
xlabel('Begin')
xticks(1:4)

ylabel('In-field firing rate (Hz)')
ax = gca;
ylim([0 ceil(ax.YLim(2))])
legend(lh, leg, 'Location', 'northeastoutside')

subplot(1,2,2)
meanData = cellfun(@mean, pkFirRate);
semData = cellfun(@semfunct, pkFirRate);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:4, meanData(g,:), semData(g,:), semData(g,:), 'Color', rgb(cols{g}));
    leg{g} = [group(g).name ' n = ' num2str(length(pkFirRate{g})) ' cells'];
end %group

xlim([0 5])
xlabel('Begin')
xticks(1:4)

ylabel('Peak firing rate (Hz)')
ax = gca;
ylim([0 ceil(ax.YLim(2))])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 2 - PLACE FIELD CENTER

figtitle = 'PlaceFieldCenter_byBeg';
figure('Name', figtitle, 'Position', [594 467 699 420])

meanData = cellfun(@mean, pfCent);
semData = cellfun(@semfunct, pfCent);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:4, meanData(g,:), semData(g,:), semData(g,:), 'Color', rgb(cols{g}));
    leg{g} = [group(g).name ' n = ' num2str(length(pfCent{g})) ' cells'];
end %group

xlim([0 5])
xlabel('Begin')
xticks(1:4)

ylabel('Place field center difference (deg)')
zero_line;
% ax = gca;
% ylim([0 ceil(ax.YLim(2))])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% FIG 3 - PK POS

figtitle = 'PlaceFieldPkPos_byBeg';
figure('Name', figtitle, 'Position', [594 467 699 420])

meanData = cellfun(@mean, pkPos);
semData = cellfun(@semfunct, pkPos);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:4, meanData(g,:), semData(g,:), semData(g,:), 'Color', rgb(cols{g}));
    leg{g} = [group(g).name ' n = ' num2str(length(pkPos{g})) ' cells'];
end %group

xlim([0 5])
xlabel('Begin')
xticks(1:4)

ylabel('Peak position difference (deg)')
zero_line;
% ax = gca;
% ylim([0 ceil(ax.YLim(2))])
legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot


%% FIG 5 - SKEWNESS

figtitle = 'PlaceFieldSkewness_byBeg';
figure('Name', figtitle, 'Position', [594 467 699 420])

meanData = cellfun(@mean, pfSkew);
semData = cellfun(@semfunct, pfSkew);

lh = nan(1,2);
leg = cell(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:4, meanData(g,:), semData(g,:), semData(g,:), 'Color', rgb(cols{g}));
    leg{g} = [group(g).name ' n = ' num2str(length(pfSkew{g})) ' cells'];
end %group

xlim([0 5])
xlabel('Begin')
xticks(1:4)

ylabel('Place field skewness')

legend(lh, leg, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot


end %function