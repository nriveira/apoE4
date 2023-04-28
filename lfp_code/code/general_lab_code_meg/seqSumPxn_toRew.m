function [sumToRew, degInFront] = seqSumPxn_toRew(pxn, rewLocs, actPos)
% function sumToRew = seqSumPxn_toRew(pxn, rewLocs, actPos)
%
% PURPOSE:
%   This function sums the posterior probability (P(x|n)) ahead of the rat
%   to the future reward location. This uses the Zheng et al. definition:
%   sum of the maximum 5 posterior probabilities from ~8 deg ahead of the
%   rat's current location to the reward lcoation.
%
% NOTE:
%   This function was written for the circle track project. It may need to
%   be modified for uses in other experiements.
%
% INPUT:
%   pxn = posterior probability distribution across the sequence event
%       (position bins x time).
%   rewLocs = reward locations on the track.
%       NOTE: stop location may be appropriate, depending on experiment.
%       Code assumes there are either 1 or 2 reward locations.
%   actPos = actual postition of the animal during sequence event.
%
% OUTPUT:
%   sumToRew = sum of max 5 posterior probabilities from ~8 deg in from of
%   	current position of the rat to the reward.
%   degInFront (optional) = number of degrees in front the cut-off starts.
%       (Code notes apprx. 8 degrees, this will give you a more exact
%       number depending on the spatial bin size used for the pxn.)
%
% MMD
% 2/2022
% Colgin Lab

spatBinSz = 360 / size(pxn,1);
degBinCtrs = spatBinSz/2:spatBinSz:360; %bin center, in deg

binsInFront = round(8/spatBinSz);
degInFront = binsInFront * spatBinSz;

if length(rewLocs) == 2
    distToRew = rad2deg(circ_dist(deg2rad(rewLocs), deg2rad(actPos))); %get dist from cur pos to next reward (in rad)
    
    nextRew = rewLocs(distToRew > 0); %find which is upcoming
else
    nextRew = rewLocs; %otherwise just use the one reward location
end %two reward locations

nextRewBin = match(nextRew, degBinCtrs); %find bin
actPosBin = match(actPos, degBinCtrs); %find bin

if nextRewBin > actPosBin+2
    pullPosBins = actPosBin+binsInFront:nextRewBin;
else
    pullPosBins = [actPosBin+binsInFront:length(degBinCtrs) 1:nextRewBin];
end %crossing 0

pullPxn = pxn(pullPosBins,:); %pull ou the appropriate bins
linPxn = pullPxn(:); %linearize
maxValsToRew = maxk(linPxn,5); %find 5 max

sumToRew = sum(maxValsToRew);

end %function