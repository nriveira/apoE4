%function enrichproj_visualizepfs(day)
% function enrichproj_visualizepfs(day)
%
% PURPOSE:
%   To visualize the place fields outputted by the "get_circtrack_pfs"
%   function.
%
% INPUT:
%   Day struct, with the circRateMap.
%
% OUTPUT:
%       Figures. Saved in specified save folder.
%
% MM Donahue
% 10/2019
% Colgin Lab

saveDir = 'E:\Rat381\pfs_figures';
cd(saveDir);
tmpallrm = NaN(4,72);

condns = {'begin1', 'begin2', 'begin3', 'begin4'};

rmBinSz = 5;
minFr = 0.25; %Hz, to be considered active in the first place
minPfLen = 18; %degrees
minPkFr = 1.5; %Hz, needs a min in-field peak firing rate of 1.5 Hz

greycol = rgb('LightGray');
for g = 1
    for r = 1
        for d = 1
            for u = 1:length(day(d).begin(1).unit)
                figtitle = ['Day' num2str(d) 'Unit' num2str(u)];
                figure('Position', [344 558 1511 420], 'Name', figtitle);

                for b = 1:4 %get all of the rms for this unit to find the max
                    tmpallrm(b,:) = day(d).begin(b).unit(u).SmCircRM;
                end %begin

                maxall = max(tmpallrm(:)); %this will allow us to scale the yaxis

                for b = 2
                    tmpallrmSM(b,:) = day(d).begin(b).unit(u).circRM;
                end

                for b = 2
                    typeInd = b;
                    subtitle = ['Begin ' num2str(b) ': ' condns{typeInd}];

                    hold on;
                    subplot(1,4,b);
                    plot(5:5:360,tmpallrm(b,:));
                    hold on;
                    plot(5:5:360,tmpallrmSM(b,:), 'k');
                    title (subtitle)
                    xlim ([0 360]);
                    xlabel ('Degrees');
                    if maxall >= minPkFr
                        ylim ([0 ceil(maxall)]);
                    else
                        ylim ([0 1.5]);

                    end
                    if b == 1
                        ylabel ('Firing rate (Hz)');
                    end

                    pf = get_circtrack_pfs(tmpallrm(b,:), rmBinSz, minFr, minPfLen);


                    if ~isempty(pf)
                        for p = 1:length(pf)
                            if pf(p).pkFr >= minPkFr
                                line([pf(p).pkPos pf(p).pkPos], [0 ceil(maxall)],  'Color','red','LineStyle','--');

                                if pf(p).radPos(1) < pf(p).radPos(end)
                                    x = [pf(p).radPos(1) pf(p).radPos(1) pf(p).radPos(end) pf(p).radPos(end)];
                                    y = [0 ceil(maxall) ceil(maxall) 0];
                                    createpatch = patch(x, y, [greycol]);
                                    alpha(0.1);
                                else %need to account for place fields that wrap around the 0-360 degree mark
                                    x = [0 0 pf(p).radPos(end) pf(p).radPos(end)];
                                    y = [0 ceil(maxall) ceil(maxall) 0];
                                    createpatch = patch(x, y, [greycol]);
                                    alpha(0.1);
                                    hold on;
                                    x = [pf(p).radPos(1) pf(p).radPos(1) 360 360];
                                    y = [0 ceil(maxall) ceil(maxall) 0];
                                    createpatch = patch(x, y, [greycol]);
                                    alpha(0.1);
                                end %if
                            end %if pk Fr high enough
                        end %for all pfs
                    end %if there are place fields

                end %begins
                saveas(gcf, figtitle, 'epsc');
                saveas(gcf, figtitle, 'fig');
            end %unit
        end %day
    end %rat
end %group
%end %fxn