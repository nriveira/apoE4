function uID = get_unit_ID(TTName)
% function uID = get_unit_ID(TTName)
% 
% PURPOSE:
%   Get the tetrode and cluster unit for a unit from the .t file name.
% 
% INPUT:
%   TTName = name for the .t file (ex: 'TT1_1.t')
%
% OUTPUT:
%  uID = [tetNum clustNum]
% 
% MMD
% 06/2022
% Colgin Lab

undInd = strfind(TTName, '_');
dotInd = strfind(TTName, '.');

tetNum = str2num(TTName(3:undInd-1));
clustNum = str2num(TTName(undInd+1:dotInd-1));

uID = [tetNum clustNum];



end %function