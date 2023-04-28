function events = detect_sequence_events_replay(pxn, spkRstr, radBinCtrs)


%% OPTIONS/INITIALIZE

sampRate = 20000; %Hz - spike sampling rate
bayesWin = 20/1000; %20 ms time window, as in Hwaun & Colgin 2019
bayesStep = 10/1000; %10 ms time step, as in Hwaun & Colgin 2019

maxJumpThr = deg2rad(35); %25 degrees - apprx. 20 cm
postSpreadThr = deg2rad(20); %12 degrees - apprx. 10 cm
timeWin = 4; %5 subsequent time bins with spike
timeStep = 1; %bin step
% distanceThr = 0.; %"distance between first and last estimated position within seqeunce was more or equal to 0.07 rad"
% jumpMeanThr = 0.5;
postSpreadThresh = 10;

%% IDENTIFY POTENTIAL EVENTS FROM PXN

onset = []; %for storing potenial event inds
offset = [];

for ii = 1:timeStep:size(pxn,2)
    range = ii:ii+timeWin-1;
    if max(range)>size(pxn,2)
        break
    end %min range fits in pxn
    
    cutPxn = pxn(:,range);
    
    emptyBin = isnan(cutPxn(1,:));
    if any(emptyBin)
        continue %to next ind
    end %any nan
    
    com = zeros(1,size(cutPxn,2)); %center of mass
    m2 = zeros(1,size(cutPxn,2));
    postSpread = zeros(1,size(cutPxn,2)); %posterior spread
    for t = 1:size(cutPxn,2)
        com(t) = circ_mean(radBinCtrs', cutPxn(:,t));
        
        for j = 1:size(pxn,1)
            m2(t) = m2(t) + circ_dist(com(t), radBinCtrs(j))^2 * cutPxn(j,t);
        end %j
        postSpread(t) = sqrt(m2(t));
    end %time bin t
    
    if mean(postSpread) > postSpreadThr
        continue
    end %too much spread
    
    jump = nan(1,length(com)-1);
    for ij = 1:length(jump)
          jump(ij) = circ_dist(com(ij+1), com(ij));
    end %get jump
    
    if max(jump) > maxJumpThr
        continue
    end %too much jump
    
    onset = [onset; ii];
    offset = [offset; range(end)];
    
end %ii - time bins

keyboard

%% COMBINE OVERLAPPING EVENTS

if size(onset,1) > 1
        isi = onset(2:end)-offset(1:end-1);
        merge = find(isi <= 0);
        
        onset(merge+1) = [];
        offset(merge) = [];
end %three are potenitally things to combine










m2 = zeros(1,size(pxn,2)); %center of mass
postSpread = zeros(1,size(pxn,2));
for t = 1:size(pxn,2)
    m2(t) = circ_var(radBinCtrs', pxn(:,t));
    %     for j = 1:size(pxn,1)
    %         m2(t) = m2(t) + circ_dist(radBinCtrs(j), com(t))^2 * pxn(j,t);
    %     end % pos bin (j)
    postSpread(t) = sqrt(m2(t));
end %time

postSpread = rad2deg(postSpread)


[S s] = circ_var(alpha, w, d, dim)

end %function