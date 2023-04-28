function [spkRstr, timeRstr] = make_spike_rasterOLD(startTm, endTm, binSz, spkTms)
% function spkRstr = make_spike_raster(startTm, endTm, binSz, spkTms)
%

spkTms = spkTms(spkTms>=startTm & spkTms<=endTm); %check
% numBins = round((endTm - startTm) / binSz);
timeRstr = startTm:binSz:endTm;
numBins = length(timeRstr);

spkRstr = zeros(1,numBins); %initialize

timePassed = spkTms - startTm;
spkInds = round(timePassed * 1/binSz);

spkInds(spkInds==0) = 1;
try
spkRstr(spkInds) = 1;
catch; keyboard; end
end %function