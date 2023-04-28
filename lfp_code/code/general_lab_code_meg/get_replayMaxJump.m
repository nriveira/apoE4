function maxJump = get_replayMaxJump(pxn, radBinCtrs)
% function maxJump = get_replayMaxJump(pxn)
%
% PURPOSE:
%   Get max jump distance between two consecutive temporal bins in the
%   replay event.
%
% INPUT:
%   pxn = deocded probability distribution across event
%   radBinCtrs = center of position bins, in radians
%
% OUTPUT:
%   maxJump = maxium jump distance between two temporal bins (in radians)
%
% MMD
% Colgin Lab
% 6/2022

com = NaN(1,size(pxn,2)); % initialize center of mass
bins2use = find(~isnan(sum(pxn,1)));
for t = 1:length(bins2use) %across time bins to use
    bin = bins2use(t); %which time bin it is
    com(bin) = circ_mean(radBinCtrs', pxn(:,bin)); %weighted circular mean based on prob distribution
end %bins - i

jumpDist = nan(1,length(com)-1);
for ji = 1:length(jumpDist)
    tmpJump = circ_dist(com(ji+1), com(ji));
    jumpDist(ji) = abs(tmpJump);
end %jumpInd

maxJump = nanmax(jumpDist);

end %function