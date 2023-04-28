function adCircTrack_x_checkRunDirection(group)
% function fmr1CircTrack_x_checkRunDirection(group)
% 
% PURPOSE:
%   For visual inspection of rad pos for the first begin of first day to 
%   check run direction.
% 
% INPUT:
%   group struct
% 
% OUTPUT:
%   F1-N: Plot of rat's radial position over time. 
%   NOTE: This function assumes the rat runs in same direction on all days
%       (and in all begins within a day, for that matter) so it just plots
%       from the first begin on the first day.
% 
% MMD
% 6/2021
% Colgin Lab

%% INITIALIZE

d = 1; %assuming same direction all days
b = 1;

%% MAKE PLOTS
for g = 1
    for r = 1:length(group(g).rat)
    radPos = group(g).rat(r).day(d).begin(b).radPos;
    
    figure
    plot(radPos((1:100:end),2))
    title([group(g).rat(r).name])
    ylim([0 360])    
    ylabel('RadPos')
    xlabel('Time')
    end %rat
end %group

end %function