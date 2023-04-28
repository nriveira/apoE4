function corrTP = replayWeighCorr(pxn, spatBinSz, varagin)
% function corrTP = replayWeighCorr(pxn, spatBinSz, bayesStep)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% DOES THIS WORK IF POS IS CIRCULAR VARIABLE I NEED TO THINK ABOUT HOW TO
% FIX IT
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


bayesStep = 10/1000; %10 ms
if nargin > 2
    bayesStep = varagin{4};
end

degBinCtrs = spatBinSz/2:spatBinSz:360; %centers of each spatial bin, in deg

iProbi = 0; %initialize
iProbP = 0;
iProbT = 0;
%get weighted means

for xi = 1:size(pxn,2)
    for yi = 1:size(pxn,1)
        iProbi = iProbi + pxn(yi,xi);
        iProbP = iProbP + pxn(yi,xi) * degBinCtrs(yi);
        iProbT = iProbT + pxn(yi,xi) * bayesStep*xi;
    end %y axis
end %x axis

mT = iProbT / iProbi;
mP = iProbP / iProbi;
%get covariance
covNum = 0; %numerator for calculating covariance
for xi = 1:size(pxn,2)
    for yi = 1:size(pxn,1)
        tmpNum = pxn(yi,xi)*((bayesStep*xi)-mT)*(degBinCtrs(yi)-mP);
        covNum = covNum + tmpNum;
    end %yi
end %xi

covTP = covNum / iProbi; %get weighted covariance for time and position
covT = cov([bayesStep * (1:size(pxn,2))]);
covP = cov(degBinCtrs);
%calculate weighted corr
corrTP = covTP/(sqrt(covT*covP));


end %function