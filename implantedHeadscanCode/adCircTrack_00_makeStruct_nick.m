function group = adCircTrack_00_makeStruct_nick
% function group = fmr1CircTrack_00_makeStruct
% 
% PURPOSE:
%   Quick function to run through all of the itemized struct-building
%   functions for the FMR1 Circle Track project.
% 
% NOTE:
%   Struct will save with data attached. See function to add a "nickname"
%   or note the the struct title.
%   See bottom of function for notes/updates on the group struct.
% 
% MMD
% 07/2021
% Colgin Lab
addpath(genpath('C:\Users\nrive\Projects\Colgin Lab\apoE4\lfp_code'));
addpath(genpath('C:\Users\nrive\Projects\Colgin Lab\GENERAL_LAB_CODE'));

curDir = pwd;
structNickName = 'nick'; 

tic
group = adCircTrack_1_buildDataStruct; %calls fmr1CircTrack_0_
group = adCircTrack_2_attachPfs(group);
group = adCircTrack_3_tagLaps_nick(group);

toc

cd(pwd)

tmpDate = clock;
cellDate = cell(1,3);
cellDate{1} = num2str(tmpDate(1));

if tmpDate(2) < 10
    cellDate{2} = ['0' num2str(tmpDate(2))];
else
    cellDate{2} = num2str(tmpDate(2));
end %if we need to add a 0 to month
if tmpDate(3) < 10
    cellDate{3} = ['0' num2str(tmpDate(3))];
else
    cellDate{3} = num2str(tmpDate(3));
end %add 0 to day
strDate = [cellDate{1} cellDate{2} cellDate{3}];

save(['dataStruct_postFxn3_' structNickName strDate], 'group')
close all
cd(curDir)

% UPDATE 7/21: Added detect sequences code
% UPDATE 8/13/21: Changed from 5 deg to 4 deg bin size
% UPDATE 9/23/21: Detect sequences/decoding with two potenital methods
% UPDATE 9/28/21: Changed min firing rate for unit to be included to 0.5 Hz
%   Added day 3 for Rat326Z
%   Fixed seq times in detect_sequence_events
% UPDATE 9/30/21: Changed min firing rate for unit bac to 1 Hz
%   Changed code to detect SWRs using Ernie method (DetectRipples_v4)
% UPDATE 10/1/21: Changed fxn 5 to only add the sleep info (spikes, coords,
%   etc). Fxn 6 now is for adding replay events via the population events
%   code and the detected SWRs from the LFP.
% UPDATE 12/3/21: Cells have been cut in sleep on 3/08 for Rat 330.
% UPDATE 12/8/21: Cells have been cut in sleep on 3/15 for Rat 330.

end %function