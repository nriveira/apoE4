function corrProj_7_4_waterFallPlots(cellRegion)
% function corrProj_7_4_waterFallPlots(cellRegion)
%
% PURPOSE:
%   Function to plot "waterfall plots for CA1 and MEC
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
binSizes = [0.005 0.010 0.050]; %ms
stateNames = {'RUN', 'REM', 'NREM','LT'}; 
lagTtls = {'+/- 5 s', '+/- 1 s'}; 
d=1;
for reg = 1%:2
    stXCorrsXState = cell(3,4); %binSize x state
    rmCcs = []; 
    for cp = 1:length(cellRegion(reg).cellPair)
        if(numel(cellRegion(reg).cellPair(cp).day) < d)
            continue;
        end
        try
            rmCcs = [rmCcs cellRegion(reg).cellPair(cp).day(d).state(1).rmCorrCoeffs];  %#ok
        catch
            rmCcs = [rmCcs cellRegion(reg).cellPair(cp).day(d+1).state(1).rmCorrCoeffs];  %#ok
        end
        for s = 1:4
            for b = 1:3 %bin
                tmpB = b+1; %wasn't plotted for the 2ms bin
                try
                    tmpStXCorr = cellRegion(reg).cellPair(cp).day(d).state(s).stXCorr{tmpB}; 
                catch
                    tmpStXCorr = cellRegion(reg).cellPair(cp).day(d+1).state(s).stXCorr{tmpB}; 
                end
                tmpStXCorr = tmpStXCorr ./ mean(tmpStXCorr); 
                stXCorrsXState{b,s} = [stXCorrsXState{b,s}; tmpStXCorr']; 
            end
        end
    end
    
   
    [~,srtInds] = sort(rmCcs);
    figure('name', regNames{reg}, 'Position', [ 280    81   780   592]);
    for b = 1:3
        for s = 1:4
            
            % Plot +/-5 s
            subplot(3,4,(b-1)*4+s);
            hold on;
            tmpXCorrs = stXCorrsXState{b,s};
            tmpXCorrs = tmpXCorrs(srtInds,:);
            xVals = 0:binSizes(b):(size(tmpXCorrs,2)-1)*binSizes(b);
            xVals = xVals - 5;
            imagesc(xVals, 1:size(tmpXCorrs,1),tmpXCorrs);
            axis xy
            colormap hot
%             caxis([.5 2])
            xlim([-1 1]);
            ylim([0 size(tmpXCorrs,1)]);
            set(gca, 'FontName', 'Arial')
            
            if b == 1
                title({stateNames{s}; lagTtls{1}})
            end
            if s == 1
                ylabel({[num2str(binSizes(b)*1000) ' ms Bins']; 'Cell Pair ID'})
            end
            if b == 3
                xlabel('Time Lag (s)'); 
            end
            
        end
        
%         % Plot +/-1 s
%         subplot(3,4,(b-1)*4+4);
%         midInd = median(1:size(tmpXCorrs,2));
%         oneSecNumInds = 1/binSizes(b);
%         numInds = 2*oneSecNumInds+1;
%         xVals = 0:binSizes(b):(numInds-1)*binSizes(b);
%         xVals = xVals - 1;
%         imagesc(xVals, 1:size(tmpXCorrs,1), tmpXCorrs(:,midInd-oneSecNumInds:midInd+oneSecNumInds))
%         xlim([-5 5]);
%         axis xy
%         colormap hot
% %         caxis([.5 2]); 
%         ylim([0 size(tmpXCorrs,1)]); 
%         set(gca, 'FontName', 'Arial')
%         if b == 1
%             title({stateNames{s}; lagTtls{2}})
%         end
%         if s == 1
%             ylabel({[num2str(binSizes(b)*1000) ' ms Bins']; 'Cell Pair ID'})
%         end
%         if b == 3
%             xlabel('Time Lag (s)');
%         end
        
    end %bin Size
   
   
    
end




end %fnctn
                
                