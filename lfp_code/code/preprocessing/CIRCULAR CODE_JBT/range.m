function y = range(x,dim)
%RANGE  Sample range.
%   Y = RANGE(X) returns the range of the values in X operating along the
%   first non-singleton dimension.
%
%   Y = RANGE(X,DIM) returns the range of the values in X operating along
%   dimension DIM.
%
%   RANGE treats NaNs as missing values, and ignores them.
%
%   See also MAX, MIN, STD.

%   Copyright 2015 The MathWorks, Inc.

if nargin < 2
   y = max(x) - min(x);
else
   y = max(x, [], dim) - min(x, [], dim);
end
