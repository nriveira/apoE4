function boutSpkTms = findSpksInBouts(spkTms, boutTms)
% function boutSpkTms = findSpksInBouts(spkTms, boutTms)
%
% INPUT:
%   spkTms = spike times for the cell
%   boutTms = n x 2 matrix with bout times

boutSpkTms = [];

for b = 1:size(boutTms,1)
    startTm = boutTms(b,1);
    endTm = boutTms(b,2);
    
    boutSpkTms = spkTms(spkTms>=startTm & spkTms<=endTm);
end %bout times

end %function