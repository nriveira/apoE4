function r2 = circLinRegress_r2(actVal, regVal)
% function r2 = circLinRegress_r2(actVal, regVal)
% 
% PURPOSE:
%   Compute the r2 value of a circular linear regression line. 
% 
% INPUT:
%   actVal = actual values that the regression line was fit to, in rad (y)
%   regVal = values of the regression line, in rad (f)
% 
% OUTPUT:
%   r2 = r2 value of the regression line
%       SSres = sum((yi - fi)^2)
%       SStot = sum((yi - mean(y))^2)
%       r2 = 1- SSres/SStot 
% 
% MD
% 3/2022
% Colgin Lab

if length(actVal) ~= length(regVal)
    error('Actual and fit values must be the same length')
end %check 

if size(actVal,1) == 1
    actVal = actVal';
end %rearrange if needed

res = zeros(1,length(actVal));
tot = zeros(1,length(actVal));

meanY = circ_mean(actVal);

for i = 1:length(actVal)
    res(i) = circ_dist(actVal(i), regVal(i))^2;
    tot(i) = circ_dist(actVal(i), meanY)^2;
end %bins

SSres = sum(res);
SStot = sum(tot);

r2 = 1 - SSres/SStot;

if r2<0
    r2 = NaN;
end



end %function