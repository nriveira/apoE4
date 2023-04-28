function propClose = get_replayPropClose(pxn, calphase, radBinCtrs, varargin)
% function propClose = get_replayPropClose(pxn, calphase, radBinCtrs, closeDef)
% 
% PURPOSE:
%   Determine the proportion of the deocded probability distribution that
%   is close to the calculated fit line.
% 
% INPUT:
%   pxn - decoded probability distribution across replay event
%   calphase = calculated fit line across replay event
%   radBinCtrs = centers of position bins, in radians
%   closeDef (optional) = definition of what is considered "close", in rad.
%       Default is 0.4 rad (apprx. 20 cm on the circle track).
% 
% OUTPUT:
%   propClose = proportion of decoded probability distribution that is 
%     close to the calculated fit line.
% 
% MMD
% Colgin Lab
% 6/2022

closeDef = 0.4; %rad - equivalent to 20 cm on the circle track
if nargin > 3
closeDef = varargin{1};
end

bins2use = find(~isnan(sum(pxn,1)));

if length(bins2use) ~= length(calphase)
    error('Calphase size must match number of decoded bins')
end

sumClose = 0; %initialize
sumAll = 0;
for t = 1:length(bins2use)
    calVal = calphase(t);
    [~,binMin] = min(abs(circ_dist(radBinCtrs, (calVal - closeDef))));
    [~,binMax] = min(abs(circ_dist(radBinCtrs, (calVal + closeDef))));
    
    if binMax > binMin
        pullBins = pxn(binMin:binMax,bins2use(t));
    else
        pullBins = [pxn(binMin:end,bins2use(t)); pxn(1:binMax,bins2use(t))];
    end %whether cross 0
    
    sumClose = sumClose + sum(pullBins);
    sumAll = sumAll + sum(pxn(:,bins2use(t)));
end %time bins

propClose = sumClose / sumAll;


end %function