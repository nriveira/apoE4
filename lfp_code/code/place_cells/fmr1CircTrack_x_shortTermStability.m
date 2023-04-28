function fmr1CircTrack_x_shortTermStability(group)
% function fmr1CircTrack_x_shortTermStability(group)
% 
% PURPOSE:
%   Examine short-term stability of place cells by computing spatial
%   correlation between first and last half of each begin. Inspired by
%   differences seen in short-term stability in WT and KO mice in Arbab,
%   Pennartz et al. 2018. The 10-min begin sessions for circle track are
%   pretty short so I'm not sure how robust or useful this will be.
%
% INPUT:
%   group data struct
%
% OUTPUT:
%   F1: Spatial correlation between the first and second halves of each
%       begin (/session) as a line plot and bar plot.
% 
% NOTE:
%   Complementary function fmr1CircTrack_x_lapRasterPlots(group) could be
%   used for illustrative purposes if there are differences in short term
%   stability between groups.
%
% MM Donahue
% 5/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 0;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\placeCellProperties';

prepForStats = 0;

cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};

%% INITIALIZE

shortTermCorr = cell(2,4); %group x session

%for get_ratemap code
spatBinSz = 4; %degrees
velFilt = 1;
durCrit = 1;

jitter = 0.1; %for dotplot

curDir = pwd;

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            for b = 1:4
                radPos = group(g).rat(r).day(d).begin(b).radPos;
                coords = group(g).rat(r).day(d).begin(b).coords; %going to use it a lot
                
                begStart = coords(1,1);
                begEnd = coords(end,1);
                
                halfInd = round(size(coords,1)/2);
                begHalf = coords(halfInd,1);
                
                for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                    
                    spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                    
                    spkTmsFirst = find(spkTms < begHalf);
                    spkTmsLast = find(spkTms >= begHalf); 
                    
                    rateMapFirst = get_ratemap_circtrack(spkTms(spkTmsFirst), coords(1:halfInd,:), radPos(1:halfInd,:), spatBinSz, velFilt, durCrit);
                    rateMapLast = get_ratemap_circtrack(spkTms(spkTmsLast), coords(halfInd:end,:), radPos(halfInd:end,:), spatBinSz, velFilt, durCrit);
                    
                    scMat = corrcoef(rateMapFirst, rateMapLast);
                    
                    shortTermCorr{g,b} = [shortTermCorr{g,b} scMat(2)];

                end %unit
            end %begin

        end %day
    end %rat
end %group
cd(saveDir)
keyboard

%% FIGS

figtitle = 'ShortTermStability';

figure('Name', figtitle, 'Position', [396 529 1158 449])

subplot(1,2,1)

tmpMean = cellfun(@nanmean, shortTermCorr);
tmpSEM = cellfun(@nansemfunct, shortTermCorr);

lh = NaN(1,2);
for g = 1:2
    hold on;
    
    lh(g) =  errorbar(1:4, tmpMean(g,:), tmpSEM(g,:), 'Color', cols{g}, 'LineWidth', 1);
     legText{g} = [groupNames{g} ': n = ' num2str(length(shortTermCorr{g,1})) ' cells'];
end %group

xlim([0 5])
xticks(1:4)
xlabel('Session')
ylabel('Spatial correlation')
legend(lh, legText, 'Location', 'southwest')
title('Spatial correlation of ratemaps from first and second half of session')

subplot(1,2,2)

tmpMean = tmpMean';
tmpSEM = tmpSEM';

bgraph = bar(1:4, tmpMean, 'FaceColor', 'Flat');
hold on;

ngroups = 4;
nbars = 2;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for g = 1:2
    bgraph(g).CData = rgb(cols{g});
    
    x = (1:ngroups) - groupwidth/2 + (2*g-1) * groupwidth / (2*nbars);
    er = errorbar(x, tmpMean(:,g), tmpSEM(:,g), '.');
    er.Color = [0 0 0];
end %group

xticks(1:4)
xlabel('Session')
ylabel('Spatial correlation')
legend(legText, 'Location', 'northeastoutside')

ax = gca;
subplot(1,2,1)
ylim([ax.YLim])

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% STATS

if prepForStats == 1
    
    statShortTermCorr = [];
    
    for g = 1:2
        for u = 1:length(shortTermCorr{g,1}) %same num of cells in all sessions for each group
            
            tmpStat = g;
            for b = 1:4
                tmpStat = [tmpStat shortTermCorr{g,b}(u)];
            end %begin
            statShortTermCorr = [statShortTermCorr; tmpStat];
            
        end %unit
    end %group
    
    keyboard
end %stats
cd(curDir)

end %function