function cpCorrVals = corrProj_7_2_1_Day1LT_scatterPlots(cellRegion)
% function corrProj_7_2_scatterPlots(cellRegion)
%
% PURPOSE:
%   Function to plot scatterplots and calculate correlations for the spike-time cross-correlation
%   (2ms bins, summed across middle 5 bins) correlated with the rate map correlation coefficients.
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
regSymbs = {'.', '^'}; 
% stateCols = {'Green', 'Gold', 'Purple'};
stateCols = {[0.0 0.5 0.0], [1.0 0.54688 0.0],[0.29297 0.0 0.50781], [0 0 0.23]};
figure('name', 'Scatterplots by Region and State', 'Position', [364 168 1051 663]); 
d=1;
for reg = 1%:2
    cpCorrVals = []; 
    for cp = 1:length(cellRegion(reg).cellPair)
        for s = 1:4
            
            if s == 1
                try
                    rmCc = cellRegion(reg).cellPair(cp).state(s).rmCorrCoeffs; 
                catch
                    continue;
                end
            end
            
            try
                midSum(s) = cellRegion(reg).cellPair(cp).state(s).midSum(1);%#ok 
            catch
                continue
            end
        end %state
        
        cpCorrVals = [cpCorrVals; rmCc midSum ]; %#ok 
        
    end %cell pair

    
    % PLOTTING
    for s = 1:4
        
        %Calculate correlation coefficients
        [rVal, pVal] = corr(cpCorrVals(:,1), cpCorrVals(:,s+1));
        fprintf('\n%s, Day 1',stateNames{s})
        var_table = table(cpCorrVals(:,1), cpCorrVals(:,s+1),'VariableNames',{'SpaceCorr','SpikeCorr'});
        mdl = fitlm(var_table,'responsevar','SpikeCorr')
        %Do the actual plotting
        subplot(2,3,(reg-1)*3+s);
        hold on; 
%         plot(cpCorrVals(:,1), cpCorrVals(:,s+1), regSymbs{reg}, 'Color', rgb(stateCols{s}));
        plot(cpCorrVals(:,1), cpCorrVals(:,s+1), '.', 'Color', (stateCols{s}),'markersize',10);
        
        lineCoefs = polyfit(cpCorrVals(:,1), cpCorrVals(:,s+1), 1);
        slope = lineCoefs(1); 
        yInt = lineCoefs(2); 
        lnCoords = [-.3*slope+yInt 1*slope+yInt]; 
        ln = line([-.3 1], lnCoords); 
        set(ln, 'Color', [0 0 0], 'LineWidth', 2); 
        
        if reg == 1
            title({stateNames{s}; ['r = ' num2str(round(rVal,4)) '; p = ' num2str(round(pVal,4))]})
        else
            title(['r = ' num2str(round(rVal,4)) '; p = ' num2str(round(pVal,4))])
        end
        if s == 1
            ylabel({regNames{reg}; 'Cross-Correlation'; 'Coefficient'})
        end
        if reg == 2 && s == 2
            xlabel('Rate Map Correlation Coefficient')
        end
        set(gca, 'FontName', 'Arial');
        ylim([0 .20]); 
        xlim([-.4 1]); 
        set(gca, 'YTick', 0:.05:.20); 
        set(gca,'xtick',-0.3:0.3:0.9)
        stats(1,reg,s) = rVal; %#ok
        stats(2,reg,s) = pVal; %#ok
    end
    
end %region
fprintf('\nNumber of cell pairs: %i\n',size(cpCorrVals,1));

% Print out stats
fprintf('STATS:\n'); 
for reg = 1%:2
    fprintf('\t%s\n', regNames{reg}); 
    for s = 1:4
        fprintf('\t\t%s:\n', stateNames{s}); 
        fprintf('\t\t\tr = %0.06g\n', stats(1,reg,s)); 
        fprintf('\t\t\tp = %0.06g\n', stats(2,reg,s)); 
    end
end
    
end %fncnt