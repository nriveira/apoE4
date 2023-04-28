function shPxn = shiftPxn(pxn, actPos, spatBinSz, vargin)
% function shPxn = shiftPxn(pxn, actPos, spatBinSz, vargin)
% 
% PURPOSE:
%   Shift pxn circuluarly to center around the rats actual position.
% 
% INPUTS:
%   pxn = output from BayesianDecoder
%   actPos = actual position of the rat in degrees
%   spatBinSz = spatial bin size in degrees used for BayesianDecoder
%   normPos = new position desiginated at rat's current position (aka what
%       degree pxn is centered around). If not input, default is 180 deg
% 
% OUTPUT:
%   shPxn = shift pxn
% 
% MMD
% 9/2021
% Colgin Lab


normPos = 180;
if nargin > 3
    normPos = varargin{1};
end

actInd = round(actPos/spatBinSz);

%shift pxn to actPos
shiftDeg = normPos - actPos;
shiftVal = round(shiftDeg/spatBinSz); % the divisor is the radial bin size in degrees

shPxn = circshift(pxn, shiftVal); %shift to be 180


end %function