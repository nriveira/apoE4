function quickNormRm(vargin)
% function quickNormRm(dayFolder)
%
% PURPOSE:
%    To plot ratemaps for each unit from one day to check if cells are
%    spatially modulated.
%
% INPUT:
%    dayFolder = directory address for folder from this day
%           ex. dayFolder = 'E:\FMR1_SWRS\RAW_DATA\Rat326Z\2021-06-11';
%           optional - can also run if currently in the correct folder
%
% OUTPUT:
%    Figures. Also saved as .pngs in day folder.
%
%  NOTE:
%    Folder must contain "TTList.txt" with accepted clusters.
%
%    All ratemaps are normalized to max firing rate in that day. See
%    function quickRawRm for ratemaps normalized within each begin window.
%
%    This function will work for big square and circle track (any
%    environment that is 1m x 1m square or eqivalent for 2d ratemap).
%
% MMD
% 06/2021
% Colgin Lab

%% INITIALIZE

if nargin == 0
    dayFolder = pwd;
end %use input or not

xBnds = [0 100]; %cm - big square is 1m x 1m, circle track has 1m diameter
yBnds = [0 100]; %cm
spatBinSz = 5; %cm
plotOrNot = 0;
velFilt = 1;
durCrit = 1;

numBins = 20; %along both x and y axis

boxSize = 9; %for smoothing

dashInd = strfind(dayFolder, '\');
dashInd = dashInd(end);
dayName = dayFolder(dashInd+1:end);

colMap = [1 1 1; jet(length(-0.01:0.01:1))]; %jet with grey nan option

cd(dayFolder)
clustlist = ReadFileList('TTList.txt');

spOpts = [1:4; 5:8; 9:12; 13:16; 17:20; 21:24]; %subplot options

allRms = nan(numBins, numBins, 4, length(clustlist));
uMaxFr = nan(1,length(clustlist)); %for adding to plot

%% GET DATA

fprintf('\tGetting ratemaps for %d units\n', length(clustlist))
for b = 1:4
    cd(['begin' num2str(b)])
    
    coords = read_in_coords('VT1.nvt', xBnds(2), yBnds(2));
    
    for u = 1:length(clustlist)
        filename = clustlist{u};
        if isfile(filename)
            spkTms = Readtfile(filename);
            spkTms = spkTms ./ 10^4; %convert to seconds
        else
            fprintf('\t\t%s - Begin %d - no spikes\n', filename, b)
            %             error('Need to write files for this tetrode!') %I forget to do this a lot
            spkTms = [];
        end
        
        
        rateMap = get_2d_ratemap(spkTms, coords, xBnds, yBnds, spatBinSz, plotOrNot, velFilt, durCrit);
        smRateMap = smooth_2d_ratemap(rateMap, boxSize);
        allRms(:,:,b,u) = smRateMap;
    end %unit
    
    cd ../
end %begin

%% MAKE FIGS

numFigs = ceil(length(clustlist)/6); %how many figure do we need - 6 units per fig

uCntr = 1;
for f = 1:numFigs
    figTitle = [dayName '_NormRateMaps_Figure' num2str(f)];
    
    figure ('Position', [680 63 832 915], 'Name', figTitle);
    
    for figUCntr = 1:6
        if uCntr <= length(clustlist)
            uName = clustlist{uCntr};
            undInd = strfind(uName, '_');
            uName = [uName(1:undInd-1) '\' uName(undInd:end)];
            
            unitRms = allRms(:,:,:,uCntr);
            maxFr = max(unitRms(:));
            
            for b = 1:4
                subplot(6,4, spOpts(figUCntr,b))
                
                plotMap = unitRms(:,:,b) ./ maxFr;
                plotMap(isnan(plotMap)) = -0.01;
                
                rmHndl = imagesc(1:numBins*spatBinSz, 1:numBins*spatBinSz, plotMap);
                axis square
                colormap(colMap)
                caxis([-0.01 1])
                set(gca, 'FontSize', 8)
                
                if b == 1
                    ylabel({[uName]; ['Max = ' num2str(round(maxFr)) ' Hz']; 'Position (cm)'});
                end %if first begin window
                
                if figUCntr == 6 || uCntr == length(clustlist) %if last unit in this fig
                    xlabel('Position (cm)');
                end
                
                if figUCntr == 1
                    title(['Begin ' num2str(b)]);
                end %if first unit
                
            end %begin
            uCntr = uCntr + 1;
        end %if there are still units to plot
    end %units in this fig
    
    saveas(gcf, figTitle, 'png')  %save png of this fig in day folder
    
end %figures


end %function