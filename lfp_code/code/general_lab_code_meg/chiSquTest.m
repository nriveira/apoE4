function [chiSquStat, df, pVal] = chiSquTest(contData, testData)
% function [chiSquStat, df, pVal] = chiSquTest(compData, obsData)
% 
% PURPOSE:
%   Perform a chi square test using the control data and test data.
% 
% INPUTS:
%   contData: control data, value you are comparing against. This input is
%       used to compute the expected proportion Ex: WT data.
%   testData: data you are comparing to the control to determine if
%       proportions are similar. Ex: KO data.
% 
% OUTPUTS:
%   chiSquStat: chi^2 statistic
%   df: degrees of freedon
%   pVal: p value
% 
% MMD
% 7/2021
% Colgin Lab


if length(contData) ~= length(testData)
   error('Inputs must be the same length!') 
end

df = length(contData) - 1; 

propAll = zeros(length(contData));
propAll = contData/sum(contData);

expData = propAll .* sum(testData);

chiSquStat = 0; %initialize
for i = 1:length(contData)
    
    tmpStat = (testData(i) - expData(i)).^2 ./ expData(i);
    chiSquStat = chiSquStat + tmpStat; %sum across
    
end

pVal = 1 - chi2cdf(chiSquStat,df);

end %function