function pp = phase_precession_entireCircTrack(spkTms, radPos, coords, lfpStruct)


%% OPTIONS

runThresh = 5; %cm/s to include spikes

%% INITIALIZE

instRs = get_runspeed(coords);
smRs = smooth_runspeed(instRs);

%get rid of spikes where rat was traveling less than threshold

badSpks = []; %initialize
for st = 1:length(spkTms)
    instSpd = smRs(find(smRs(:,1)<=spkTms(st), 1, 'Last'),2);
    if instSpd < runThresh
        badSpks = [badSpks st];
    end %if too slow
    
end %spkTms
  
spkTms(badSpks) = []; %new spike times with slow spikes removed

pp.spkTms = spkTms;
pp.spkPhis = zeros(length(spkTms),1);
pp.spkPos = zeros(length(spkTms),1);
pp.spkAngPos = zeros(length(spkTms),1);
pp.stats = zeros(2,1);

angToDegScale = (pi*100)/360; %because track has 100cm diameter;

phiTms = get_theta_phase_times(lfpStruct); %get theta phase times (pk, fa, tr, ri)
phiVctr = get_asym_theta_phi_vector(lfpStruct); %get a theta phase vector, interpolated from phase times

for st = 1:length(spkTms)
    if spkTms(st) > lfpStruct.ts(1)
        spkLfpInd = find(lfpStruct.ts<=spkTms(st), 1, 'Last'); %get the closest LFP ind
        
        pp.spkPhis(st) = phiVctr(spkLfpInd); %get the theta phase at which spike occurred
        
        posInd = find(radPos(:,1)<=spkTms(st), 1, 'Last'); %closest position index
        pp.spkAngPos(st) = radPos(posInd,2); %radial position at which spike occurred
        pp.spkPos(st) = angToDegScale .* radPos(posInd,2);
    end
end %spkTms

%I need to do circular-circular regression here
[r2 pVal] = circ_corrcc(pp.spkPos,pp.spkPhis); %calculate circular/linear regression of spike phases/positions
pp.stats = [r2 pVal];


end %function