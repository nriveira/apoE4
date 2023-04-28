function group = adCircTrack_3_tagLaps_jayanth(group)
% 
% PURPOSE
%  Functions add indices to main data structure that point to times which
%  the rat was completing a lap based on the crossing od radpos from 0 to 360 degree
%  (hence one lap). This code does not consider any reward indices
%  information to get half laps.
%
% INPUT:
%  group = main data structure for the project
%
% OUTPUT:
%  A subfield within each 'begin' sub-structure that contains the indices for the start and stop of full laps
%
% Jayanth 
% 2/10/23
% Colgin Lab

% plotCheck = 1; %Way to doublecheck that paths are covering the expected portions of the track
%              Makes 2 figures for each begin window (for each day/rat/group).
%savePlots = 1;
saveDir = 'E:\Colgin Lab\resultsApril2023_AD_WT\lapPlots';
cd(saveDir);

for g = 1:2
    
    if g == 1

        fprintf('Group %d\n', g);

        for r = 1:length(group(g).rat)
            %for r = 1
            fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);

            for d = 1:length(group(g).rat(r).day)
                %for d = 3
                fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));

                for b = 1:4
                    fprintf('\t\t\tBegin %d\n', b);

                    radPos = group(g).rat(r).day(d).begin(b).radPos;
                    coords = group(g).rat(r).day(d).begin(b).coords;

                    radPosUn = unwrap(deg2rad(radPos(:,2)));
                    radPosUn = radPosUn + abs(min(radPosUn));
                    numLaps = abs(round(round(radPosUn(1) - radPosUn(end))/(2*pi)));
                    LapindTmp = [];

                    for i = 1:numLaps

                        LapindTmp(i) = max(find(radPosUn > (radPosUn(1) - 2*pi*i)));

                        if i == 1
                            group(g).rat(r).day(d).begin(b).lapInds(i,:) = [1 LapindTmp(i)];
                            group(g).rat(r).day(d).begin(b).lapTms(i,:) = [radPos(1,1) radPos(LapindTmp(i),1)];
                            figure;
                            plot(coords(1:LapindTmp(i),2), coords(1:LapindTmp(i),3));
                            figtitle = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_begin' num2str(b) '_lap' num2str(i)];
                            saveas(gcf, figtitle, 'epsc');
                            saveas(gcf, figtitle, 'fig');
                            saveas(gcf, figtitle, 'png');
                            close;
                            continue;
                        end

                        if i == numLaps
                            LapindTmp(i) = length(radPosUn);
                            group(g).rat(r).day(d).begin(b).lapInds(i,:) = [LapindTmp(i-1)+1 LapindTmp(i)];
                            group(g).rat(r).day(d).begin(b).lapTms(i,:) = [radPos(LapindTmp(i-1)+1,1) radPos(LapindTmp(i),1)];
                            figure;
                            plot(coords(LapindTmp(i-1)+1:LapindTmp(i),2), coords(LapindTmp(i-1)+1:LapindTmp(i),3));
                            figtitle = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_begin' num2str(b) '_lap' num2str(i)];
                            saveas(gcf, figtitle, 'epsc');
                            saveas(gcf, figtitle, 'fig');
                            saveas(gcf, figtitle, 'png');
                            close;
                            continue;
                        end

                        group(g).rat(r).day(d).begin(b).lapInds(i,:) = [LapindTmp(i-1)+1 LapindTmp(i)];
                        group(g).rat(r).day(d).begin(b).lapTms(i,:) = [radPos(LapindTmp(i-1)+1,1) radPos(LapindTmp(i),1)];
                        figure;
                        plot(coords(LapindTmp(i-1)+1:LapindTmp(i),2), coords(LapindTmp(i-1)+1:LapindTmp(i),3));
                        figtitle = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_begin' num2str(b) '_lap' num2str(i)];
                        saveas(gcf, figtitle, 'epsc');
                        saveas(gcf, figtitle, 'fig');
                        saveas(gcf, figtitle, 'png');
                        close;

                    end % numLaps

                end %begin

            end %day

        end %rat
    
    end % for clockwise (my ad rats)

    if g == 2
        fprintf('Group %d\n', g);

        for r = 1:length(group(g).rat)
            %for r = 1
            fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);

            for d = 1:length(group(g).rat(r).day)
                %for d = 3
                fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));

                for b = 1:4
                    fprintf('\t\t\tBegin %d\n', b);

                    radPos = group(g).rat(r).day(d).begin(b).radPos;
                    coords = group(g).rat(r).day(d).begin(b).coords;

                    radPosUn = unwrap(deg2rad(radPos(:,2)));
                    radPosUn = radPosUn + abs(min(radPosUn));
                    numLaps = abs(round(round(radPosUn(1) - radPosUn(end))/(2*pi)));
                    LapindTmp = [];

                    for i = 1:numLaps

                        LapindTmp(i) = max(find(radPosUn < (radPosUn(1) + 2*pi*i)));

                        if i == 1
                            group(g).rat(r).day(d).begin(b).lapInds(i,:) = [1 LapindTmp(i)];
                            group(g).rat(r).day(d).begin(b).lapTms(i,:) = [radPos(1,1) radPos(LapindTmp(i),1)];
                            figure;
                            plot(coords(1:LapindTmp(i),2), coords(1:LapindTmp(i),3));
                            figtitle = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_begin' num2str(b) '_lap' num2str(i)];
                            saveas(gcf, figtitle, 'epsc');
                            saveas(gcf, figtitle, 'fig');
                            saveas(gcf, figtitle, 'png');
                            close;
                            continue;
                        end

                        if i == numLaps
                            LapindTmp(i) = length(radPosUn);
                            group(g).rat(r).day(d).begin(b).lapInds(i,:) = [LapindTmp(i-1)+1 LapindTmp(i)];
                            group(g).rat(r).day(d).begin(b).lapTms(i,:) = [radPos(LapindTmp(i-1)+1,1) radPos(LapindTmp(i),1)];
                            figure;
                            plot(coords(LapindTmp(i-1)+1:LapindTmp(i),2), coords(LapindTmp(i-1)+1:LapindTmp(i),3));
                            figtitle = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_begin' num2str(b) '_lap' num2str(i)];
                            saveas(gcf, figtitle, 'epsc');
                            saveas(gcf, figtitle, 'fig');
                            saveas(gcf, figtitle, 'png');
                            close;
                            continue;
                        end

                        group(g).rat(r).day(d).begin(b).lapInds(i,:) = [LapindTmp(i-1)+1 LapindTmp(i)];
                        group(g).rat(r).day(d).begin(b).lapTms(i,:) = [radPos(LapindTmp(i-1)+1,1) radPos(LapindTmp(i),1)];
                        figure;
                        plot(coords(LapindTmp(i-1)+1:LapindTmp(i),2), coords(LapindTmp(i-1)+1:LapindTmp(i),3));
                        figtitle = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_begin' num2str(b) '_lap' num2str(i)];
                        saveas(gcf, figtitle, 'epsc');
                        saveas(gcf, figtitle, 'fig');
                        saveas(gcf, figtitle, 'png');
                        close;

                    end % numLaps

                end %begin

            end %day

        end %rat
    
    end % for counterclockwise (meg's WT rats)

end %group

end % function

% close all;
% P = {};
% for i = 1:length(LapindTmp)
%     if i == 1
%         figure;
%         plot(coords(i:LapindTmp(i), 2), coords(i:LapindTmp(i), 3)); 
%         radPosTmp{i,1} = radPos(1:LapindTmp(i),2);
%         x = 1:length(radPosTmp{i,1});
%         P{i,1} = polyfit(x,radPosTmp{i,1},1);
%         figure;
%         scatter(x , radPosTmp{i,1},'b','*'); 
%         yfit = P{i,1}(1)*x+P{i,1}(2);  % P(1) is the slope and P(2) is the intercept
%         hold on;
%         plot(x,yfit,'r-.')
%         continue;
% 
%     elseif i == length(LapindTmp)
%         figure;
%         plot(coords(LapindTmp(i): length(radPos), 2), coords(LapindTmp(i): length(radPos), 3)); 
%         radPosTmp{i,1} = radPos(LapindTmp(i):length(radPos),2);
%         x = 1:length(radPosTmp{i,1});
%         P{i,1} = polyfit(1:length(radPosTmp{i,1}),radPosTmp{i,1},1);
%         figure;
%         scatter(x , radPosTmp{i,1},'b','*'); 
%         yfit = P{i,1}(1)*x+P{i,1}(2);  % P(1) is the slope and P(2) is the intercept
%         hold on;
%         plot(x,yfit,'r-.')
%         continue;
%     end
% 
%     figure;
%     plot(coords(LapindTmp(i-1)+1:LapindTmp(i), 2), coords(LapindTmp(i-1)+1:LapindTmp(i), 3));
%     radPosTmp{i,1} = radPos(LapindTmp(i-1)+1:LapindTmp(i),2);
%     x = 1:length(radPosTmp{i,1});
%     P{i,1} = polyfit(1:length(radPosTmp{i,1}),radPosTmp{i,1},1);
%     figure;
%     scatter(x , radPosTmp{i,1},'b','*'); 
%     yfit = P{i,1}(1)*x+P{i,1}(2);  % P(1) is the slope and P(2) is the intercept
%     hold on;
%     plot(x,yfit,'r-.')
% end
