function adCircTrack_x_checkRewLoc(group)


for g = 1
    for r = 1
        for d = 1:length(group(g).rat(r).day)
            figtitle = [group(g).name '_' group(g).rat(r).name '_D' num2str(d)];
            figure('Name', figtitle)
            rewLocs = group(g).rat(r).day(d).rewLocs;
            rewLocs(rewLocs==0) = 360;
            hold on;
            
            for b = 1:4
                 radPos = group(g).rat(r).day(d).begin(b).radPos;
                plot(radPos(:,2))
            end %begin
            for rl = 1:length(rewLocs)
                line([0 20000], [rewLocs(rl) rewLocs(rl)], 'Color', 'Black', 'LineStyle', '--', 'LineWidth', 1.5) 
                end %rew loc
        end %day
    end %rat
end %group


end %function