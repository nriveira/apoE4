function [r2, pVal, slope, int] = circ_lin_regression(linData, circData) 
% Circ data in radians!
% This is a variation of CZ's code, but that code shuffles and shifts pxn
keyboard



[para,~,p] = circ_lin_regress(linData, circData); %CZ code

calphase = 2*pi*para(1,1)*linData+para(1,2);
slope = para(1,1)*2*pi;
%     calphase = calphase*360/2/pi;
over = calphase>=2*pi;
under = calphase<0;
calphase(over) = calphase(over)-2*pi;
calphase(under) = calphase(under)+2*pi; %from -2pi to 2pi
keyboard
%compute residual
yRes = zeros(1,length(linData));
yTot = zeros(1,length(linData));
for i = 1:length(linData)
    yFit = calphase(i);
    yRes(1,i) = circdistance(circData(i), yFit, 1);
    yTot(1,i) = circdistance(circData(i), circ_mean(circData), 1);
end

SSresid = sum(yresidual);
SStotal = sum(ytotal);

r2 = 1 - SSresid/SStotal;
if r2<0
    r2 = NaN;
end




end %function