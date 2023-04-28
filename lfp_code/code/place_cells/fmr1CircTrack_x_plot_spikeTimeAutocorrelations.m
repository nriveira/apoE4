function fmr1CircTrack_x_plot_spikeTimeAutocorrelations(group)
% function fmr1CircTrack_x_plot_spikeTimeAutocorrelations(group)
%
% PURPOSE:
%   Plot spike time autocorrelations for cells from WT and KO rats.
%
% INPUTS:
%   group struct
%
% OUTPUTS:
%   Fn-n (optional): Spike time autocorrelations across all units in day.
%   Fn-n (optional): Spike time autocorrelations across all units in rat.
%   F1: Spike time autocorrelation across all units for both genotypes.
%
% OPTIONS:
%   See code for options to plot by individual rats and days. Options to
%   save figs.
%
% MMD
% 07/2021
% Colgin Lab

%% OPTIONS

plotByDay = 0; %1 to plot mean for every day, 0 to not
plotByRat = 0; %1 to plot mean for every rat, 0 to not

saveOrNot = 1;
saveDir = 'E:\resultsFeb2023_AD_WT\spikeTimeAutocorrelations';
curDir = pwd;

cols = {'Blue', 'Red'};

%% INITIALIZE

timeLag = 1; %s
binSz = 10/1000; %10 ms

autoCorr = cell(2,1);

minFr = 1; %Hz, to be included

cd(saveDir)

%% GET DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\t%s\n', group(g).rat(r).name)
        autoCorrByRat = [];
        
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            dayEdges = cell(1,4);
            for b = 1:4
                timeStart = group(g).rat(r).day(d).begin(b).coords(1,1);
                timeEnd = group(g).rat(r).day(d).begin(b).coords(end,1);
                
                tmpEdges = timeStart:binSz:timeEnd;
                dayEdges{b}  = tmpEdges;
            end %get bin edges
            
            autoCorrByDay = [];
            
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                fprintf('\t\t\tUnit %d/%d\n', u, length(group(g).rat(r).day(d).xAllBeginUnitInfo))
                if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap) < minFr
                    continue
                end %if it doesn't reach firing rate threshold
                autoCorrForUnit = [];
                
                for b = 1:4
                    spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                    if isempty(spkTms)
                        continue %to next
                    end %no spk tms
                    tmpEdges = dayEdges{b};
                    if length(tmpEdges) > 2^16 %check - hist counts will not work if there are too many bins
                        keyboard
                    end
                    
                    tmpBinSpks = histcounts(spkTms, tmpEdges);
                    tmpAutoCorr = xcorr(tmpBinSpks, timeLag/binSz, 'coeff');
                    autoCorrForUnit = [autoCorrForUnit; tmpAutoCorr];
                end %begin
                
                autoCorrByRat = [autoCorrByRat; mean(autoCorrForUnit,1)];
                autoCorrByDay = [autoCorrByDay; mean(autoCorrForUnit,1)];
                autoCorr{g} = [autoCorr{g}; mean(autoCorrForUnit,1)];
            end %unit
            
            %Make within day figure
            if plotByDay == 1
                if isempty(autoCorrByDay)
                    continue %to next day
                end %if there is not data
                figName = [group(g).name '_' group(g).rat(r).name '_' group(g).rat(r).day(d).name];
                figure('Name', figName)
                
                plotData = mean(autoCorrByDay,1);
                plotData(ceil(length(plotData)/2)) = NaN;
                
                plot(linspace(-timeLag, timeLag, length(plotData)), plotData, 'Color', rgb(cols{g}), 'LineWidth', 1.5)
                
                xlim([-timeLag timeLag])
                xticks(-timeLag:1:timeLag)
                xlabel('Time (s)')
                
                ttl = {[group(g).name ' - ' group(g).rat(r).name ' - ' group(g).rat(r).day(d).name]; ['n = ' num2str(size(autoCorrByDay,1)) ' cells']};
                title(ttl)
                
                ylim([0 .2])
                ylabel('Correlation coefficient')
                y_zero_line;
                
                if saveOrNot == 1
                    cd([saveDir '\byDay'])
                    saveas(gcf, figName, 'epsc')
                    saveas(gcf, figName, 'png')
                    saveas(gcf, figName, 'fig')
                end %save option
            end %if plot by day
            
        end %day
        
        % Make all days for this rat figure
        if plotByRat == 1
            figName = [group(g).name '_' group(g).rat(r).name '_allDays'];
            figure('Name', figName)
            
            plotData = mean(autoCorrByRat,1);
            plotData(ceil(length(plotData)/2)) = NaN;
            plot(linspace(-timeLag, timeLag, length(plotData)), plotData, 'Color', rgb(cols{g}), 'LineWidth', 1.5)
            
            xlim([-timeLag timeLag])
            xticks(-timeLag:1:timeLag)
            xlabel('Time (s)')
            
            ttl = {[group(g).name ' - ' group(g).rat(r).name]; ['n = ' num2str(size(autoCorrByRat,1)) ' cells']};
            title(ttl)
            
            ylim([0 .2])
            ylabel('Correlation coefficient')
            y_zero_line;
            
            if saveOrNot == 1
                cd([saveDir '\byRat'])
                saveas(gcf, figName, 'epsc')
                saveas(gcf, figName, 'png')
                saveas(gcf, figName, 'fig')
            end %save option
        end %plot by rat
        
    end %rat
end %group

cd(saveDir)
keyboard

%% FIG 1(ish) 

figtitle = 'SpikeTimeAutoCorrelations';
figure('Name', figtitle, 'Position', [576 528 861 420])

for g = 1:2
    subplot(1,2,g)
    
    plotData = mean(autoCorr{g},1);
    plotData(ceil(length(plotData)/2)) = NaN;
    plot(linspace(-timeLag, timeLag, length(plotData)), plotData, 'Color', rgb(cols{g}), 'LineWidth', 1.5)
    
    xticks(-timeLag:1:timeLag)
    xlabel('Time (s)')
    
    ttl = {[group(g).name]; ['n = ' num2str(size(autoCorr{g},1)) ' cells']};
    title(ttl)
    
    ylabel('Correlation coefficient')
    y_zero_line;
end %group

same_axes

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

cd(curDir)

end %function