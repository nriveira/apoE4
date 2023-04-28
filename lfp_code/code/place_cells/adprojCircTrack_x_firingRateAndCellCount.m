function adprojCircTrack_x_firingRateAndCellCount(group)
% function fmr1CircTrack_x_firingRateAndCellCount(group)
%
% PURPOSE:
%  To plot firing rate histograms for each group/day/rat with unit counts
%  in the subplot titles.
%
% INPUT:
%  group = project uber data struct
%
% OUTPUT:
%   F1-n: Figure for each rat with subplots for individual days showing
%       firing rate histograms.
%
% JB Trimper
% 01/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1;
saveDir = 'E:\resultsFeb2023_AD_WT\placecell_props\firingRates';

%% INITIALIZE

curDir = pwd;
cd(saveDir)

%% GET DATA/MAKE PLOTS

for g = 1:2
    
    for r = 1:length(group(g).rat)
        figtitle = group(g).rat(r).name;
        figure('Position', [205.8 183.4 859.2 594.4], 'Name', figtitle);
        
        numU = zeros(1,length(group(g).rat(r).day));
        
        for d = 1:length(group(g).rat(r).day)
            numU(d) = length(group(g).rat(r).day(d).xAllBeginUnitInfo);
            spkCnts = zeros(1, length(group(g).rat(r).day(d).xAllBeginUnitInfo));
            
            sesDur = 0;
            for b = 1:4
                startBegin = group(g).rat(r).day(d).begin(b).coords(1,1);
                endBegin = group(g).rat(r).day(d).begin(b).coords(end,1);
                sesDur = sesDur + endBegin - startBegin;
                
                for u = 1:length(group(g).rat(r).day(d).begin(b).unit)
                    spkCnts(u) = spkCnts(u) + length(group(g).rat(r).day(d).begin(b).unit(u).spkTms);
                end
                
            end %begin
            
            firingRates = spkCnts / sesDur;
            
            subplot(1,5,d);
            histogram(firingRates, 'BinWidth', 0.25);
            title(['Day' num2str(d) ' (n = ' num2str(numU(d)) ')']);
            if mod(d,5) == 1
                ylabel('Unit Count');
            end
            if d>length(group(g).rat(r).day)-5
                xlabel('Firing Rate (Hz)');
            end
            yBnds = get(gca, 'Ylim');
            yMax = yBnds(2)+.2*yBnds(2);
            ylim([0 yMax]);
            xBnds = get(gca, 'XLim');
            if xBnds(2) <= 2
                xBnds(2) = 2;
            end
            xlim([-0.05 xBnds(2)]);
            
        end %day
        
        if saveOrNot == 1
            saveas(gcf, figtitle, 'epsc');
            saveas(gcf, figtitle, 'fig');
            saveas(gcf, figtitle, 'png');
        end %save or not
        
    end %rat
end %group

cd(curDir)
end %fnctn