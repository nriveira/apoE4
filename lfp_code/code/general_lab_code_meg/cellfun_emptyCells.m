function resultArray = cellfun_emptyCells(funcHandle, cellArray, custFill)
% function resultArray = cellfun_emptyCells(funcHandle, cellArray, custFill)
% 
% PURPOSE:
%   This function will complete the cellfun function across a cell array
%   that has empty cells. It can be used when you are looking to get a
%   single value from applying the function to each cell. For example, when
%   you are trying to get the standard error of the mean across all cells. 
%   It's marginally quicker than writing this out each time.
% 
% INPUT:
%   cellArray = cell array containing the data
%   funcHandle = function handle to be applied to every cell (include @)
%       ex. @semfunct
%   custFill (optional) = what you want to fill into the empty cells. If
%       left blank, function will put Nan.
% 
% OUTPUT:
%   resultArray = matrix of results of applying funcHandle across all
%       cells in cellArray, with custFill used where there is no data.
% 
% MMD
% 01/2021
% Colgin Lab

if nargin < 3
    custFill = NaN;
end %check custom fill
try
tmpResult = cellfun(funcHandle, cellArray, 'UniformOutput', false);
catch; keyboard; end
isEmpt = cellfun(@isempty, tmpResult);

tmpResult(isEmpt) = {custFill};
resultArray = cell2mat(tmpResult);

end %function