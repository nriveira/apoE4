function pp = phase_precession_circtrack_pf(spkTms, radPos, coords, pf, lfpStruct)
%
% PURPOSE:
%   Function assesses phase precession, providing spike phases by position and phase-precession stats
%    NOTE: The direction that the rat was running (i.e., clockwise vs. counterclockwise [180 toward 360 vs
%          360 toward 180) will impact the results here. You may have to adjust lines 53-57 (see notes next to those lines)
%
% INPUT:
%     spkTms = vector of spike times for the unit
%     coords = matrix of rat's coordinates by time (output of 'read_in_coords')
%         pf = structure for each place-field (output of 'get_circtrack_pfs')
%  lfpStruct = lfp structure with narrow theta filterested lfps
%           John originally used asym_theta_phi_vector, but I changed it to
%           just used angles of hilbert transformed singals.
%
% OUTPUT:
%    pp = phase precession (length(pp) = # pfs for this unit) struct with subfields:
%          1) pp.passSpkTms (spk times for each pass)
%          2) pp.passSpkPhis (spk theta phases for each pass in degrees)
%          3) pp.passSpkPos (spk positions in cm)
%          4) pp.passSpkPosRel (spk positions in cm, relative to pf start, for each pass)
%          5) pp.passSpkAngPos (angular spk positions)
%          6) pp.passSpkAngPosRel (angular spk positions, relative to pf start, for each pass)
%          7) pp.passNormSpkPos (within-field spike positions normalized between 0 and 1)
%          8) pp.passPPStats (phiPrecess stats; matrix of slope, intercept, R^2, and pVal - USES RELATIVE SPK POSITION)
%
% OPTIONS:
%   This function only includes spike times where the rat was moving
%   >5cm/s. Run threshold can be changed internally.
%
%
% JB Trimper - Edited by MMD for as through pf, not by pass
% 09/2019 - Edited 06/2021
% Colgin Lab

%% OPTIONS

runThresh = 5; %cm/s to include spikes - 10 in Feng, Silva, Foster 2015

bound = 2;

th1 = 6; %lower bound for theta
th2 = 10; %upper bound for theta

%% INITIALIZE

angToDegScale = (pi*100)/360; %because track has 100cm diameter;

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

% phiVctr = angle(hilbert(-lfpStruct.narrowThetaLfp))*180/pi+180;%theta peaks = 0 and 360 - from CZ code
% phiTms = get_theta_phase_times(lfpStruct); %get theta phase times (pk, fa, tr, ri)
% phiVctr = get_asym_theta_phi_vector(lfpStruct); %get a theta phase vector, interpolated from phase times


bpth = fftbandpass(lfpStruct.data, lfpStruct.Fs, th1-2,th1,th2,th2+2);%theta
phiVctr = angle(hilbert(-bpth))*180/pi+180;%theta peaks = 0 and 360

for p = 1:length(pf) %for each place-field
    
    pp(p).spkTms = [];
    pp(p).spkPhis = [];
    pp(p).spkPos = []; %actual track position, in cm
    pp(p).spkPosRel = []; %position in cm relative to pf start
    pp(p).spkAngPos = []; %actual rad pos
    pp(p).spkAngPosRel = []; %rad pos relative to pf start
    pp(p).normSpkPos = []; %positions in pf normalized between 0 and 1
    pp(p).stats = zeros(4,1);
    
    %If rat was running from 360 toward 180, use this line
    %             startPos = pf(p).radPos(end);
    
    %If rat was running from 180 toward 360, use this line
    startField = pf(p).radPos(1);
    endField = pf(p).radPos(end);
    
    fieldSz = rad2deg(circ_dist(deg2rad(endField), deg2rad(startField)));
    
    if endField > startField
        crossZero = 0;
        tmpfieldSz = endField - startField;
    else
        crossZero = 1;
        tmpfieldSz = (360 - startField) + endField;
    end
    if round(fieldSz) ~= tmpfieldSz
        keyboard
    end %check
    
    for st = 1:length(spkTms) %for each spike
        
        inField = 0; %initialize as not in this pf
        
        %         spkPos = radPos(find(radPos(:,1) <= spkTms(st), 1, 'Last'),2); %get rad pos when spike occured
        spkPos = radPos(match(spkTms(st), radPos(:,1)),2);
        if isempty(spkPos)
            keyboard
            %             spkPos = radPos(find(radPos(:,1) >= spkTms(st), 1, 'First'),2);
        end
        
        if crossZero == 0
            if spkPos >= startField && spkPos <= endField
                inField = 1;
                pp(p).spkTms = [pp(p).spkTms spkTms(st)];
            end
        elseif spkPos >= 0 && spkPos <= endField || spkPos >= startField
            inField = 1;
            pp(p).spkTms = [pp(p).spkTms spkTms(st)];
            
        end %does pf cross zero
        
        
        if inField == 1
            %             spkLfpInd = find(lfpStruct.ts<=spkTms(st), 1, 'Last'); %get the closest LFP ind
            spkLfpInd = match(spkTms(st), lfpStruct.ts); %get the closest LFP ind
            %             if isempty(spkLfpInd)
            %                 spkLfpInd = find(lfpStruct.ts>=spkTms(st), 1, 'First');
            %             end
            pp(p).spkPhis = [pp(p).spkPhis phiVctr(spkLfpInd)]; %get the theta phase at which spike occurred
            
            %             posInd = find(radPos(:,1)<=spkTms(st), 1, 'Last'); %closest position index
            pp(p).spkAngPos = [pp(p).spkAngPos spkPos]; %radial position at which spike occurred
            pp(p).spkPos = [pp(p).spkPos spkPos .*angToDegScale];
            
            tmpRelAng =  abs(rad2deg(circ_dist(deg2rad(spkPos), deg2rad(startField))));
            tmpSpkPosRel = tmpRelAng .* angToDegScale; %convert angular distances of spike positions from peak to centimeter distances
            
            pp(p).spkAngPosRel = [pp(p).spkAngPosRel tmpRelAng];
            pp(p).spkPosRel = [pp(p).spkPosRel tmpSpkPosRel];
            
            tmpNormPos = tmpRelAng / fieldSz;
            pp(p).normSpkPos = [pp(p).normSpkPos tmpNormPos]; 
        end %it is in the pf
    end %spkTms
    
    %     para = circ_lin_regress(pp(p).normSpkPos, deg2rad(pp(p).spkPhis), bound); %using CZ code and methods
        beta = CircularRegression(pp(p).normSpkPos, deg2rad(pp(p).spkPhis));
    %     if para(1,1) == bound || para(1,1) == - bound
    %         keyboard
    %     end
    
    %     calphase = 2*pi*para(1,1)*(0:1) + para(1,2);
    calphase = beta(1)*[0, 1] + beta(2);
    
    tmpCal = beta(1)*pp(p).normSpkPos + beta(2);
    
    r2 = circLinRegress_r2(deg2rad(pp(p).spkPhis), tmpCal);
    
    %     hold on; plot(0:1, calphase)
    
    %     slope = calphase(2) - calphase(1); %because norm pos - 0-1 (would be divide by 1)
    slope = beta(1);
    if sum(calphase < 0) < 0 && sum(calphase>=2*pi) > 0
        keyboard
    end
    
    if calphase(1) < 0 %phase off between 0 and 2*pi
        calphase = calphase + 2*pi;
    elseif calphase(1) >= 2*pi
        calphase = calphase -2*pi;
    end
    
    phaseOff = calphase(1);
    if phaseOff < 0 || phaseOff >= 2*pi
        keyboard
    end
    
    if ~isempty(pp(p).spkPhis)
        %         [R2 pval] = circ_corrcl(deg2rad(pp(p).spkPhis), pp(p).normSpkPos);
        %         pp(p).stats = [slope phaseOff R2 pval];
        pp(p).stats = [slope phaseOff r2];
    else
        pp(p).stats = [];
    end
    
    
end %placefields




end %fnctn