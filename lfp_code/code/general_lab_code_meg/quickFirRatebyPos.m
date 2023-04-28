function quickFirRatebyPos(varagin)
% function quickFirRatebyPos(dayFolder)
%

if nargin == 0
    dayFolder = pwd;
end %use input or not


spatBinSz = 4; %degrees
numBins = 360/spatBinSz; %360/5 degree bins
velFilt = 1;
durCrit = 1;

gWinDur = 30; %degrees
gWinStd = gWinDur/2; %degrees

dashInd = strfind(dayFolder, '\');
dashInd = dashInd(end);
dayName = dayFolder(dashInd+1:end);

cd(dayFolder)
clustlist = ReadFileList('TTList.txt');

% spOpts = [1:4; 5:8; 9:12; 13:16; 17:20; 21:24]; %subplot options
spOpts = [1:5; 6:10; 11:15; 16:20; 21:25; 26:30]; %subplot options

allRms = nan(numBins, 4, length(clustlist));

for u = 1:length(clustlist)
    fprintf('\tUnit %d/%d\n', u, length(clustlist))
    for b = 1:4
        
        cd(['begin' num2str(b)])
        filename = clustlist{u};
        if ~isfile(filename)
            fprintf('\t\t%s - Begin %d\n', filename, b)
            error('Need to write files for this tetrode!') %I forget to do this a lot
        end %check that files were written
        
        spkTms = Readtfile(filename);
        spkTms = spkTms ./ 10^4; %convert to seconds
        
        [radPos, coords] = read_in_circtrack_coords('VT1.nvt');
        
        rateMap = get_ratemap_circtrack(spkTms, coords, radPos, spatBinSz, velFilt, durCrit);
        smRateMap = smooth_circtrack_ratemap(rateMap, gWinDur, gWinStd);
        
        allRms(:,b,u) = smRateMap;
        
        cd ../
    end %begin
end %unit


numFigs = ceil(length(clustlist)/9); %how many figure do we need - 9 units per fig

uCntr = 1;

for f = 1:numFigs
    figTitle = [dayName '_NormRateMaps_Figure' num2str(f)];
    
    figure ('Position', [680 63 832 915], 'Name', figTitle);
    
    for figUCntr = 1:9
        if uCntr <= length(clustlist)
              uName = clustlist{uCntr};
            undInd = strfind(uName, '_');
            uName = [uName(1:undInd-1) '\' uName(undInd:end)];
            
            subplot(3,3,figUCntr)
            
            for b = 1:4
                hold on;
                plot(spatBinSz/2:spatBinSz:360, allRms(:,b,uCntr))
                
            end %begin
            uCntr = uCntr + 1;
        end %there are still units
    end %fig unit counter
    
end %figures



end %function