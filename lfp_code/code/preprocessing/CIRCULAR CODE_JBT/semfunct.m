function sem = semfunct(x,dim)
% function sem = semfunct(x,dim)
%
% SUMMARY: 
%   Function to calculate the standard error of the sample mean along the dimension specified.
%
% INPUT:
%    x = vector or matrix of values
%  dim = optional input specifying the dimension along which the standard error should be calculated
%         If left out, dim is assigned as the first non-singleton dimension. 
%         
% JBT 8/15
% Manns Lab
% Emory Univ


if nargin == 1
    allDims = size(x);
    dim = find(allDims > 1, 1, 'First');
end

sem = nanstd(x,0,dim) ./ sqrt(size(x,dim));

end %function
