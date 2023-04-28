function corrProj_7_3_linegraphs(cellRegion)
% function corrProj_7_3_linegraphs(cellRegion)
%
% PURPOSE:
%   Function to plot line graphs for each region, with a line for each state, showing relative spatial
%   phase/distance plotted against the spike-time cross correlation coefficient, summed across the middle bins.
%
% INPUT:
%   cellRegion = the output of corrProj_6... which has data for each cell pair by state and region
%
% OUTPUT:
%   Figures (scatterplots) as well as printed stats.
%
% JBT 8/2017
% Colgin Lab


regNames = {'MEC', 'CA1'};
stateNames = {'RUN', 'REM', 'NREM','LT'};
binSizes = {'+/-5ms', '+/-50ms'};
stateCols = {'Green', 'Gold', 'Purple', 'Navy'};
% stateCols = {[0.0 0.5 0.0], [1.0 0.54688 0.0],[0.29297 0.0 0.50781], [0 0 0.23]};

scaled5BinFig = figure('name', 'Normalized, 5 Bins', 'Position', [364 168 1051 663]);
unscaled5BinFig = figure('name', 'Non-normalized, 5 Bins', 'Position', [364 168 1051 663]);

scaled6BinFig = figure('name', 'Normalized, 6 Bins', 'Position', [364 168 1051 663]);
unscaled6BinFig = figure('name', 'Non-normalized, 6 Bins', 'Position', [364 168 1051 663]);

%MEC
binEdges{1,1} = 0:.25:1;
binEdges{1,2} = 0:.2:1;
%CA1
binEdges{2,1} = 0:36:180;
binEdges{2,2} = 0:30:180; %6 bins
d=2;
for reg = 1%:2
    
    fiveBinVals = cell(2,5);
    sixBinVals = cell(2,6);
    
    for cp = 1:length(cellRegion(reg).cellPair)
        if(numel(cellRegion(reg).cellPair(cp).day) < d)
            continue;
        end
        for s = 1:4
            
            if s == 1
                if reg == 1
                    try
                        relSpatDist = cellRegion(reg).cellPair(cp).day(d).state(s).relSpatPhiMag;
                    catch
                        relSpatDist = cellRegion(reg).cellPair(cp).day(d+1).state(s).relSpatPhiMag;
                    end
                else
                    relSpatDist = cellRegion(reg).cellPair(cp).day(d).state(s).relSpatDist;
                end
            end
            try
                midSum5(s) = cellRegion(reg).cellPair(cp).day(d).state(s).midSum(1);%#ok
                midSum50(s) = cellRegion(reg).cellPair(cp).day(d).state(s).midSum(2);%#ok
            catch
                midSum5(s) = cellRegion(reg).cellPair(cp).day(d+1).state(s).midSum(1);%#ok
                midSum50(s) = cellRegion(reg).cellPair(cp).day(d+1).state(s).midSum(2);%#ok
            end
        end %state
        
        %for 5 bins
        binInd = find(binEdges{reg,1}<=relSpatDist, 1, 'Last');
        fiveBinVals{1,binInd} = [fiveBinVals{1,binInd}; midSum5];
        fiveBinVals{2,binInd} = [fiveBinVals{2,binInd}; midSum50];
        
        %for 6 bins
        binInd = find(binEdges{reg,2}<=relSpatDist, 1, 'Last');
        sixBinVals{1,binInd} = [sixBinVals{1,binInd}; midSum5];
        sixBinVals{2,binInd} = [sixBinVals{2,binInd}; midSum50];
        
    end %cell pair
    
    
    
    
    %% PLOT FIVE BIN VERSION
    
    %   First, plot the non-normalized version
    figure(unscaled5BinFig);
    AVG = zeros(5,4);
    SEM = zeros(5,4);
    for bs = 1:2 %bin size (5 or 50)
        for bn = 1:5 %bin # (1-5)
            try
                AVG(bn,:) = mean(fiveBinVals{bs,bn},1);
            catch
                keyboard
            end
            try
            
                SEM(bn,:) = std(fiveBinVals{bs,bn},0,1)./sqrt(size(fiveBinVals{bs,bn},1)); %semfunct(fiveBinVals{bs,bn},1);
            catch
                keyboard
            end
        end
        
        subplot(2,2,(reg-1)*2+bs);
        hold on;
        for s = 1:4
            tmpLn = plot(1:5, AVG(:,s));
            set(tmpLn, 'Color', rgb(stateCols{s}));
            errorbar(1:5, AVG(:,s), SEM(:,s), 'Color', rgb(stateCols{s}), 'LineStyle', 'None')
        end
        set(gca, 'XLim', [.9 5.1], 'XTick', 1:5, 'FontName', 'Arial');
        if bs == 1
            ylabel({regNames{reg}; 'Cross-Correlation Coefficient'});
        end
        if reg == 1
            xlabel('Relative Spatial Phase');
            set(gca, 'XTickLabel', {'0.00 - 0.25', '0.25 - 0.50', '0.50 - 0.75', '0.75 - 1.00', '>1.00'})
            title(binSizes{bs});
        else
            xlabel('Relative Angular Distance');
            set(gca, 'XTickLabel', {'0-36', '36-72', '72-108', '108-144', '144-180'})
        end
        yBnds = get(gca, 'YLim');
        ylim([0 yBnds(2)]);
        
        maxAVG(bs,:) = max(AVG); %for scaled figure
    end
    
    %  Then, plot the normalized version
    figure(scaled5BinFig);
    for bs = 1:2 %bin size (5 or 50)
        for bn = 1:5 %bin # (1-5)
            
            for s = 1:4
                %normalize all individual values to the max average across bins for each state
                tmpNormVals = fiveBinVals{bs,bn}(:,s) ./ maxAVG(bs,s);
                
                AVG(bn,s) = mean(tmpNormVals,1);
                SEM(bn,s) =  std(tmpNormVals,0,1)./sqrt(size(tmpNormVals,1));%semfunct(tmpNormVals,1);
            end
        end
        
        subplot(2,2,(reg-1)*2+bs);
        hold on;
        for s = 1:4
            tmpLn = plot(1:5, AVG(:,s));
            set(tmpLn, 'Color', rgb(stateCols{s}));
            legend({'Open Field','','NREM','','REM','','Linear Track'})
            errorbar(1:5, AVG(:,s), SEM(:,s), 'Color', rgb(stateCols{s}), 'LineStyle', 'None')
        end
        set(gca, 'XLim', [.9 5.1], 'XTick', 1:5, 'FontName', 'Arial');
        if bs == 1
            ylabel({regNames{reg}; 'Cross-Correlation Coefficient'});
        end
        if reg == 1
            xlabel('Relative Spatial Phase');
            set(gca, 'XTickLabel', {'0.00 - 0.25', '0.25 - 0.50', '0.50 - 0.75', '0.75 - 1.00', '>1.00'})
            title(binSizes{bs});
        else
            xlabel('Relative Angular Distance');
            set(gca, 'XTickLabel', {'0-36', '36-72', '72-108', '108-144', '144-180'})
        end
        yBnds = get(gca, 'YLim');
        ylim([0 yBnds(2)]);
    end
    
    
    
    
    
    
    
    
    
    
    
    %% PLOT 6 BIN VERSION
    
    %  First, plot the non-normalized version
    figure(unscaled6BinFig);
    AVG = zeros(6,4);
    SEM = zeros(6,4);
    for bs = 1:2 %bin size (5 or 50)
        for bn = 1:6 %bin # (1-5)
            AVG(bn,:) = mean(sixBinVals{bs,bn},1);
            SEM(bn,:) = std(sixBinVals{bs,bn},0,1)./sqrt(size(sixBinVals{bs,bn},1)); %semfunct(sixBinVals{bs,bn},1);
        end
        
        subplot(2,2,(reg-1)*2+bs);
        hold on;
        for s = 1:4
            tmpLn = plot(1:6, AVG(:,s));
            set(tmpLn, 'Color', rgb(stateCols{s}));
            errorbar(1:6, AVG(:,s), SEM(:,s), 'Color', rgb(stateCols{s}), 'LineStyle', 'None')
        end
        set(gca, 'XLim', [.9 6.1], 'XTick', 1:6, 'FontName', 'Arial');
        if bs == 1
            ylabel({regNames{reg}; 'Cross-Correlation Coefficient'});
        end
        if reg == 1
            xlabel('Relative Spatial Phase');
            set(gca,'XTickLabel', {'0.00-0.20', '0.20-0.40', '0.40-0.60', '0.60-0.80', '0.80-1.00', '>1.00'});
            title(binSizes{bs});
        else
            xlabel('Relative Angular Distance');
            set(gca,'XTickLabel', {'0-30', '30-60', '60-90', '90-120', '120-150', '150-180'});
        end
        yBnds = get(gca, 'YLim');
        ylim([0 yBnds(2)]);
        
         maxAVG(bs,:) = max(AVG); %for scaled figure
    end
    
    
    %  Then, plot the normalized version
    figure(scaled6BinFig);
    for bs = 1:2 %bin size (5 or 50)
        for bn = 1:6 %bin # (1-5)
            
            for s = 1:4
                %normalize all individual values to the max average across bins for each state
                tmpNormVals = sixBinVals{bs,bn}(:,s) ./ maxAVG(bs,s);
                
                AVG(bn,s) = mean(tmpNormVals,1);
                SEM(bn,s) = std(tmpNormVals,0,1)./sqrt(size(tmpNormVals,1)); %semfunct(tmpNormVals,1);
            end
        end
        
        subplot(2,2,(reg-1)*2+bs);
        hold on;
        for s = 1:4
            tmpLn = plot(1:6, AVG(:,s));
            set(tmpLn, 'Color', rgb(stateCols{s}));
            legend({'Open Field','','NREM','','REM','','Linear Track'})
            errorbar(1:6, AVG(:,s), SEM(:,s), 'Color', rgb(stateCols{s}), 'LineStyle', 'None')
        end
        set(gca, 'XLim', [.9 6.1], 'XTick', 1:6, 'FontName', 'Arial');
        if bs == 1
            ylabel({regNames{reg}; 'Cross-Correlation Coefficient'});
        end
        if reg == 1
            xlabel('Relative Spatial Phase');
            set(gca,'XTickLabel', {'0.00-0.20', '0.20-0.40', '0.40-0.60', '0.60-0.80', '0.80-1.00', '>1.00'});
            title(binSizes{bs});
        else
            xlabel('Relative Angular Distance');
            set(gca,'XTickLabel', {'0-30', '30-60', '60-90', '90-120', '120-150', '150-180'});
        end
        yBnds = get(gca, 'YLim');
        ylim([0 yBnds(2)]);
    end
    
    
    
end %region


end%fnctn