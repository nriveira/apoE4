Runfd_novel = {'C:\HippoReplayproject\Data\Rat_100\circulartrack\2016-04-12_16_combined';...
    'C:\HippoReplayproject\Data\Rat_112\circulartrack\2016-10-14_combine';...
    'C:\HippoReplayproject\Data\Rat_113\circulartrack\2017-01-30_combine';...
    'C:\HippoReplayproject\Data\Rat_114\circulartrack\2016-08-27_combine';...
    'C:\HippoReplayproject\Data\Rat_100\circulartrack\2016-04-13_08_combined';...
    'C:\HippoReplayproject\Data\Rat_112\circulartrack\2016-10-15_combine';...
    'C:\HippoReplayproject\Data\Rat_113\circulartrack\2017-01-31_combine';...
    'C:\HippoReplayproject\Data\Rat_114\circulartrack\2016-08-28_combine'};
trialfd = {'begin1';'begin2';'begin3';'begin4'};
for rr = 1:length(Runfd_novel);
    outdir = strcat(Runfd_novel{rr},'\ScanningEvent');
    if ~isdir(outdir)
        mkdir(outdir)
    end
    
    for ii = 1:length(trialfd)
        posfile = strcat(Runfd_novel{rr},'\',trialfd{ii},'\VT1.nvt');
        [post,posx,posy,circle] = LoadPos(posfile);
        [onset,offset, ind] = DetectScanningEvent(post,posx,posy,circle);
        
        ScanningEvent.(trialfd{ii}).onset = onset;
        ScanningEvent.(trialfd{ii}).offset = offset;
        ScanningEvent.(trialfd{ii}).ind = ind;
    end
    save(strcat(outdir,'\ScanningEvent.mat'),'ScanningEvent')
    clear ScanningEvent
end