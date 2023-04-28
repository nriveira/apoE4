function evComp = replay_decodeStructure(pxn)

maxJump = 1.4; %rad
minJump = 0.14; %rad, distance considered a "jump" or step
minDist = 0.07; %rad - for whole seq

bayesStep = 10/1000; %default as 10 ms
minBin = 0.03 / bayesStep;

spatBinSz = 2*pi / size(pxn,1);
radBinCtrs = spatBinSz/2:spatBinSz:2*pi;

timeAx = 0:bayesStep:bayesStep*size(pxn, 2);
% bins2use = find(~isnan(sum(pxn,1)));

allCom = NaN(1,size(pxn,2)); %initialize
for t = 1:size(pxn,2)
    allCom(t) = wrapTo2Pi(circ_mean(radBinCtrs', pxn(:,t))); %weighted circular mean based on prob distribution
end %time bins

nanBins = find(isnan(allCom));
for nb = 1:length(nanBins)
    allCom(nanBins(nb)) = wrapTo2Pi(circ_mean([allCom(nanBins(nb)-1) allCom(nanBins(nb)+1)]')); %interpolate the missing
end

allJump = nan(1,length(allCom)-1);
for ij = 1:length(allJump)
    allJump(ij) = circ_dist(allCom(ij+1), allCom(ij));
end %get jump

forBins = [];
revBins = [];
statBins = [];

for ii = 1:size(pxn,2)-(minBin-1)
    range = ii:ii+minBin-1;

    com = allCom(range);
    jump = allJump(ii:ii+minBin-2);

    %     if abs(max(jump)) > maxJump
    %         continue
    %     end %too big jump

    jump(jump < minJump & jump > 0) = 0; %if less than min, not considered a "real jump"
    jump(jump > -minJump & jump < 0) = 0;

    %     jci = find(diff(sign(jump)));
    jsi = sign(jump);
    jsi = jsi(jsi~=0); %ignore hovers in same spot

    if isempty(jsi)
        statBins = [statBins range];
        continue
    end %no jumps in any direction

    rDist = circ_dist(com(end), com(1));
    if rDist == 0 || isempty(jsi)
          statBins = [statBins range];
          continue
    end

    if all(jsi == sign(rDist))
        if sign(rDist) == 1
            forBins = [forBins range];
        else
            revBins = [revBins range];
        end
    end %if signs match

end %ii - time bins

allForBins = sort(unique(forBins));
allRevBins = sort(unique(revBins));
allStatBins = sort(unique(statBins));

keepStat = allStatBins(~ismember(allStatBins, union(allRevBins, allForBins)));

compCntr = 0;

tmpBins = zeros(1, size(pxn,2));
tmpBins(allForBins) = 1;
cc = bwconncomp(tmpBins);
if cc.NumObjects ~=0
    for c = 1:cc.NumObjects
        compCntr = compCntr + 1;
        evComp(compCntr).name = 'forward';
        evComp(compCntr).type = 1; %just means forward

        useInds = cc.PixelIdxList{c};
        evComp(compCntr).tInds = useInds; %time bin inds

        maxJump = max(abs(allJump(useInds(1):useInds(end)-1)));

        uwCom = unwrap(allCom(useInds)); %unwrap
        fit = polyfit(timeAx(useInds), uwCom, 1); %get fit variables for this portion
        R = corrcoef(timeAx(useInds), uwCom);
        r2 = R(2)^2;

        calphase = timeAx(useInds)*fit(1) + fit(2);

        segDist = circ_dist(uwCom(end), uwCom(1)); %distance in radians

        evComp(compCntr).r2 = r2;
        evComp(compCntr).com = allCom(useInds);
        evComp(compCntr).timeAx = timeAx(useInds);
        evComp(compCntr).calphase = wrapTo2Pi(calphase);
        evComp(compCntr).slope = fit(1);
        evComp(compCntr).dist = segDist;
        evComp(compCntr).maxJump = maxJump;
    end %con comp
end %not empty

tmpBins = zeros(1, size(pxn,2));
tmpBins(allRevBins) = 1;
cc = bwconncomp(tmpBins);
if cc.NumObjects ~=0
    for c = 1:cc.NumObjects
        compCntr = compCntr + 1;
        evComp(compCntr).name = 'reverse';
        evComp(compCntr).type = 2; %just means reverse

        useInds = cc.PixelIdxList{c};
        evComp(compCntr).tInds = useInds; %time bin inds

        maxJump = -max(abs(allJump(useInds(1):useInds(end)-1)));

        %         [r2, calphase, ~, ~, slope] = Cir_reg(pxn, radBinCtrs', timeAx, useInds);

        uwCom = unwrap(allCom(useInds)); %unwrap
        fit = polyfit(timeAx(useInds), uwCom, 1); %get fit variables for this portion
        R = corrcoef(timeAx(useInds), uwCom);
        r2 = R(2)^2;

        calphase = timeAx(useInds)*fit(1) + fit(2);

        segDist = circ_dist(uwCom(end), uwCom(1)); %distance in radians

        evComp(compCntr).r2 = r2;
        evComp(compCntr).com = allCom(useInds);
        evComp(compCntr).timeAx = timeAx(useInds);
        evComp(compCntr).calphase = wrapTo2Pi(calphase);
        evComp(compCntr).slope = fit(1);
        evComp(compCntr).dist = segDist;
            evComp(compCntr).maxJump = maxJump;
    end %con comp
end %not empty

tmpBins = zeros(1, size(pxn,2));
tmpBins(keepStat) = 1;
cc = bwconncomp(tmpBins);
if cc.NumObjects ~=0
    for c = 1:cc.NumObjects
        compCntr = compCntr + 1;
        evComp(compCntr).name = 'stationary';
        evComp(compCntr).type = 3; %just means stationary

          useInds = cc.PixelIdxList{c};
        evComp(compCntr).tInds = useInds; %time bin inds

         uwCom = unwrap(allCom(useInds)); %unwrap
%         fit = polyfit(timeAx(useInds), uwCom, 1); %get fit variables for this portion
%         R = corrcoef(timeAx(useInds), uwCom);
%         r2 = R(2)^2;

%         calphase = timeAx(useInds)*fit(1) + fit(2);

        evComp(compCntr).r2 = [];
        evComp(compCntr).com = allCom(useInds);
        evComp(compCntr).timeAx = timeAx(useInds);
        evComp(compCntr).calphase = [];
        evComp(compCntr).slope = [];
%         [r2, calphase, ~, ~, slope] = Cir_reg(pxn, radBinCtrs', timeAx, useInds');
% 
%         evComp(compCntr).r2 = r2;
%         evComp(compCntr).calphase = calphase;
%         evComp(compCntr).slope = slope;
    end %con comp
end %not empty

end %function