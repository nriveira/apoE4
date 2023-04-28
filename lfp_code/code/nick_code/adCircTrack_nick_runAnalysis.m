% Initialize variables/paths
curDir = pwd;
structNickName = 'nick_apo';
addpath(genpath("C:\Users\nrive\Projects\Colgin Lab\GENERAL_LAB_CODE"))
addpath(genpath("C:\Users\nrive\Projects\Colgin Lab\GENERAL_CODE_MEG_WROTE"))
addpath(genpath("C:\Users\nrive\Projects\Colgin Lab\apoE4\lfp_code\code"))
figure(1); clf; figure(2); clf;

%% Run functions to build group struct
group = adCircTrack_1_buildDataStruct;
group = adCircTrack_nick_getLapSpeed(group);

%% Run to save group data struct
cd('E:\Colgin Lab\data_structs')

tmpDate = clock;
if tmpDate(2) < 10
    strDate = [num2str(tmpDate(1)) '0' num2str(tmpDate(2)) num2str(tmpDate(3))];
else
    strDate = [num2str(tmpDate(1)) num2str(tmpDate(2)) num2str(tmpDate(3))];
end %if we need to add a 0 to month

save(['dataStruct_' structNickName strDate], 'group')
close all
cd(curDir)