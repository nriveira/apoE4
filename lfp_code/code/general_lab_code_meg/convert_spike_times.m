function conSpkTms = convert_spike_times(spkTms)
% function conSpkTms = convert_spike_times(spkTms)
%
% PURPOSE:
%   Convert spike times from Neuralynx sampling times to seconds. Just
%   making this a function because I always forget the conversion rate and
%   constantly have to double check it.
% 
% conSpkTms = spkTms ./ 10^4;
%
% MMD
% 06/2022
% Colgin Lab

conSpkTms = spkTms ./ 10^4;

end %function