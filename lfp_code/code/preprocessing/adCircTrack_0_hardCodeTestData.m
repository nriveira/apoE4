function group = adCircTrack_0_hardCodeTestData(group)
% function group = fmr1CircTrack_0_hardCodeTestData(group)
%
% Function adds the names, dates, reward locs, and theta tet to use to the
% group structure. It's just separated to this function vs it's parent because
% it takes less space.
%
% JBT
% 03/21
% Colgin Lab




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%                 apoE4 RATS                 %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RAT 381
group(1).rat(1).name = 'rat381';

group(1).rat(1).day(1).name = '2022-10-23'; 
group(1).rat(1).day(1).rewLocs = [225 45]; % check on this again, not sure about this value here
group(1).rat(1).day(1).thetaTet = 3; %7 units, highest amp/cleanest theta
% -->>>>> Theta tetrodes chosen by manual inspection of LFPs across all tets for each day
%          - Evaluated for # of units, clarity of theta, and amplitude of theta

group(1).rat(1).day(2).name = '2022-10-24'; %Day 2
group(1).rat(1).day(2).rewLocs = [315 135]; % E/W
group(1).rat(1).day(2).thetaTet = 3; %highest amp/cleanest theta

group(1).rat(1).day(3).name = '2022-10-25'; %Day 3
group(1).rat(1).day(3).rewLocs = [45 225]; % NW/SE
group(1).rat(1).day(3).thetaTet = 2; %highest amp/cleanest theta

group(1).rat(1).day(4).name = '2022-10-26'; %Day 4
group(1).rat(1).day(4).rewLocs = [135 315]; % NE/SW
group(1).rat(1).day(4).thetaTet = 4; % highest amp/cleanest theta

group(1).rat(1).day(5).name = '2022-11-3'; %Day 5
group(1).rat(1).day(5).rewLocs = [135 315]; % NNW/SSE
group(1).rat(1).day(5).thetaTet = 3; % highest amp/cleanest theta

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%                 WT RATS                 %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RAT 326Z

group(2).rat(1).name = 'rat326Z';

group(2).rat(1).day(1).name = '2021-06-14A'; %Day 1
group(2).rat(1).day(1).rewLocs = [0 180]; % E/W
group(2).rat(1).day(1).thetaTet = [6]; %10 units

group(2).rat(1).day(2).name = '2021-06-14B'; %Day 2
group(2).rat(1).day(2).rewLocs = [45 225]; %NE/SW
group(2).rat(1).day(2).thetaTet = [5]; %4 units, but least amount of noise

group(2).rat(1).day(3).name = '2021-06-15'; %Day 3
group(2).rat(1).day(3).rewLocs = [135 315];
group(2).rat(1).day(3).thetaTet = [5]; %EXPLAIN HERE


%%
% 
% group(1).rat(1).day(10).name = '2020-11-23'; %Day 10
% group(1).rat(1).day(10).rewLocs = [157.5 337.5]; % WNW/ESE
% group(1).rat(1).day(10).thetaTet = 2; %only option
% 
% group(1).rat(1).day(11).name = '2020-11-27'; %Day 11
% group(1).rat(1).day(11).rewLocs = [90 270]; %N/S
% group(1).rat(1).day(11).thetaTet = 2; %2 units, cleanest (others not good looking)
% 
% group(1).rat(1).day(12).name = '2020-11-28'; %Day 12
% group(1).rat(1).day(12).rewLocs = [0 180]; %E/W
% group(1).rat(1).day(12).thetaTet = 2; %3 tetNums, cleanest
% 
% group(1).rat(1).day(13).name = '2020-11-30'; %Day 13
% group(1).rat(1).day(13).rewLocs = [45 225]; %NE/SW
% group(1).rat(1).day(13).thetaTet = 4; %5 units, cleanest, high amp
% 
% %% RAT 330
% group(1).rat(2).name = 'rat330';
% group(1).rat(2).day(1).name = '2021-03-04'; %Day 1
% group(1).rat(2).day(1).rewLocs = 45; % NE (only 1 rew loc)
% group(1).rat(2).day(1).thetaTet = 3; %8 units, tied for most units, cleanest theta
% 
% group(1).rat(2).day(2).name = '2021-03-08'; %Day 2 FOR MEG for Emma this is different
% group(1).rat(2).day(2).rewLocs = [135 315]; % NW/SE
% group(1).rat(2).day(2).thetaTet = 3; %selected by Emma
% 
% group(1).rat(2).day(3).name = '2021-03-15'; %Day 3 FOR MEG for Emma this is different
% group(1).rat(2).day(3).rewLocs = [135 315]; % NW/SE
% group(1).rat(2).day(3).thetaTet = 4; %selected by Emma
% 
% 
% 
% 
% 
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                         %
% %                  WT RATS                %
% %                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% %% RAT 113
% % group(1).rat(1).name = 'rat113';
% % group(1).rat(1).day(1).name = '2017-01-29'; %Day 1
% % group(1).rat(1).day(1).rewLocs = [0 180]; % E/W
% % group(1).rat(1).day(1).thetaTet = 12; %4 units (least #), but by far the cleanest theta
% 
% %% RAT 326Z
% 
% group(1).rat(1).name = 'rat326Z';
% group(1).rat(1).day(1).name = '2021-06-14A'; %Day 1
% group(1).rat(1).day(1).rewLocs = [0 180]; % E/W
% group(1).rat(1).day(1).thetaTet = [6]; %10 units
% 
% group(1).rat(1).day(2).name = '2021-06-14B'; %Day 2
% group(1).rat(1).day(2).rewLocs = [45 225]; %NE/SW
% group(1).rat(1).day(2).thetaTet = [5]; %4 units, but least amount of noise
% 
% group(1).rat(1).day(3).name = '2021-06-15'; %Day 3
% group(1).rat(1).day(3).rewLocs = [135 315];
% group(1).rat(1).day(3).thetaTet = [5]; %EXPLAIN HERE









