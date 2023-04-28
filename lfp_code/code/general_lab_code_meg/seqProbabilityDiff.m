function probDiff = seqProbabilityDiff(shiftPxn, spatBinSz, varagin)
% function probDiff = seqProbabilityDiff(shiftPxn, spatBinSz, trackDiameter)
% 
% PURPOSE:
%   Calculate the probability difference as defined by Feng, Silva, &
%       Foster 2015. Calculated from the decoded probable location of the
%       animal across time in the sequence event. Positive differences
%       imply theta seqeunces sweeping in running direction of the animal.
%       See reference for more detailed information.
% 
% INPUTS:
%   shiftPxn = pxn during seqeunce event shifted so rat's actual position
%       is at 180 deg (see shiftPxn function)
%   spatBinSz = spatial bin size in deg for pxn
%   trackDiameter (optional) = default is 100 cm (standard circular track
%       for Colgin Lab. Include if different track size is used.
% 
% OUTPUT:
%   probDiff = computed difference between the pxn in quadrants (I and III)
%       and (II and IV). See reference paper for more information.
% 
% MMD
% 9/2021
% Colgin Lab

trackDiameter = 100;
if nargin > 2
    trackDiameter = varagin{3};
end

degCmConv = trackDiameter*pi / 360;
degBinCtrs = spatBinSz/2:spatBinSz:360;

cutY = [size(shiftPxn,1)/2 - ceil(degCmConv*50)/spatBinSz:(size(shiftPxn,1)/2 + ceil(degCmConv*50)/spatBinSz-1)]; %take deocded probability within 50 cm of actual position

cutYPxn = shiftPxn(cutY,:);

if rem(size(cutYPxn,2),2) == 0 %even
    q1 = cutYPxn(1:size(cutYPxn,1)/2, size(cutYPxn,2)/2+1:end);
    q2 = cutYPxn(1:size(cutYPxn,1)/2, 1:size(cutYPxn,2)/2);
    q3 = cutYPxn(size(cutYPxn,1)/2+1:end, 1:size(cutYPxn,2)/2);
    q4 = cutYPxn(size(cutYPxn,1)/2+1:end, size(cutYPxn,2)/2+1:end);
else %odd
    midCol = cutYPxn(:,ceil(size(cutYPxn,2)/2));
    q1 = cutYPxn(1:size(cutYPxn,1)/2, ceil(size(cutYPxn,2)/2)+1:end) + midCol(1:size(cutYPxn,1)/2)/2;
    q2 = cutYPxn(1:size(cutYPxn,1)/2, 1:floor((size(cutYPxn,2)/2))) + midCol(1:size(cutYPxn,1)/2)/2;
    q3 = cutYPxn(size(cutYPxn,1)/2+1:end, 1:floor((size(cutYPxn,2)/2))) + midCol(size(cutYPxn,1)/2+1:end)/2;
    q4 = cutYPxn(size(cutYPxn,1)/2+1:end, ceil(size(cutYPxn,2)/2)+1:end) + midCol(size(cutYPxn,1)/2+1:end)/2;
end %even or odd number of  time bins

probDiff = ((sum(q2(:))+sum(q4(:))) - (sum(q1(:))+sum(q3(:))))/sum([sum(q1(:)) sum(q2(:)) sum(q3(:)) sum(q4(:))]);

end %function