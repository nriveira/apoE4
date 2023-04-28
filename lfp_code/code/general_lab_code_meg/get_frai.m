function frai = get_frai(spkTms)
% function frai = get_frai(spkTms)
% 
% PURPOSE:
%   Get the firing rate asymmetry index from spikes that occured during a
%   pass through a place field. See Mehta et al. 2000 for reference.
% 
% INPUT:
%   spkTms = spike times through a place field pass.
% 
% OUTPUT:
%   frai = firing rate asymmetry index.
% 
% MMD
% 2/2022
% Colgin Lab

halfSpks = ceil(length(spkTms)/2);
FR1 = halfSpks/(spkTms(halfSpks) - spkTms(1));
FR2 = halfSpks/(spkTms(end) - spkTms(halfSpks));
frai = (FR1-FR2)/(FR1+FR2);

end %function