function corrTP = seqWeighCorr(shiftPxn, spatBinSz, varagin)
% function corrTP = seqWeighCorr(shiftPxn, spatBinSz, trackDiameter, bayesStep)
%
% PURPOSE:
%   To calculate the theta seqeunce weighted correlation between time (T)
%       and position (P) during seqeunce events, as in Feng, Silva, &
%       Foster 2015.
%
% INPUT:
%   shiftPxn = pxn during seqeunce event shifted so rat's actual position
%       is normalized to 180 deg (see shiftPxn function)
%   spatBinSz = spatial bin size in deg for pxn
%   trackDiameter (optional) = default is 100 cm (standard circular track
%       for Colgin Lab. Include if different track size is used.
%   bayesStep (optional) = time step for the Bayesian decoder in seconds. 
%       If not included, default is 10 ms
%
% OUTPUT:
%   corrTP = weighted correlation between time and position across sequence
%       event
%
% MMD
% 9/2021
% ColginLab

trackDiameter = 100;
if nargin > 2
    trackDiameter = varagin{3};
end

bayesStep = 10/1000; %10 ms
if nargin > 3
    bayesStep = varagin{4};
end

degCmConv = trackDiameter*pi / 360;
degBinCtrs = spatBinSz/2:spatBinSz:360; %centers of each spatial bin, in deg

cutY = [size(shiftPxn,1)/2 - ceil(degCmConv*50)/spatBinSz:(size(shiftPxn,1)/2 + ceil(degCmConv*50)/spatBinSz-1)]; %take decodec probability within 50 cm of actual position

iProbi = 0; %initialize
iProbP = 0;
iProbT = 0;
%get weighted means
for xi = 1:size(shiftPxn,2)
    for yi = cutY
        iProbi = iProbi + shiftPxn(yi,xi);
        iProbP = iProbP + shiftPxn(yi,xi) * degBinCtrs(yi);
        iProbT = iProbT + shiftPxn(yi,xi) * bayesStep*xi;
    end %y axis
end %x axis

mT = iProbT / iProbi;
mP = iProbP / iProbi;
%get covariance
covNum = 0; %numerator for calculating covariance
for xi = 1:size(shiftPxn,2)
    for yi = cutY
        tmpNum = shiftPxn(yi,xi)*((bayesStep*xi)-mT)*(degBinCtrs(yi)-mP);
        covNum = covNum + tmpNum;
    end %yi
end %xi

covTP = covNum / iProbi; %get weighted covariance for time and position
covT = cov([bayesStep * (1:size(shiftPxn,2))]);
covP = cov(degBinCtrs(cutY));
%calculate weighted corr
corrTP = covTP/(sqrt(covT*covP));

end %function