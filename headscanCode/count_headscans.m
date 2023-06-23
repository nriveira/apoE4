function output = count_headscans(headscan_struct)
    apoE3 = headscan_struct(contains({headscan_struct.group},'apoE3'));
    apoE4 = headscan_struct(contains({headscan_struct.group},'apoE4'));
    figure(1); clf; figure(2); clf

    for s = 1:3
        apoE3_session = apoE3(strcmp([apoE3.session],num2str(s)));
        apoE4_session = apoE4(strcmp([apoE4.session],num2str(s)));

        apoE3_headscans = zeros(length(apoE3_session),1);
        apoE4_headscans = zeros(length(apoE4_session),1);
        apoE3_headscans_dur = zeros(length(apoE3_session),1);
        apoE4_headscans_dur = zeros(length(apoE4_session),1);

        for h = 1:length(apoE3_session)
            apoE3_headscans(h) = apoE3_session(h).headscan.num_headscans;
            apoE3_headscans_dur(h) = diff(apoE3_session(h).headscan.filt_pss(h,:)) / 30;
        end

        for h = 1:length(apoE4_session)
            apoE4_headscans(h) = apoE4_session(h).headscan.num_headscans;
            apoE4_headscans_dur(h) = diff(apoE4_session(h).headscan.filt_pss(h,:)) / 30;
        end
 
        figure(1);
        subplot(3,2,2*(s-1)+1);
        histogram(apoE3_headscans,10);
        title(['apoE3 S' num2str(s)])
        xlabel('Number of headscans')
        ylabel('Count')

        subplot(3,2,2*s)
        histogram(apoE4_headscans,10);
        title(['apoE4 S' num2str(s)])
        xlabel('Number of headscans')
        ylabel('Count')

        figure(2);
        subplot(3,2,2*(s-1)+1);
        histogram(apoE3_headscans_dur,10);
        title(['apoE3 S' num2str(s)])
        xlabel('Headscan Duration [s]')
        ylabel('Count')

        subplot(3,2,2*s)
        histogram(apoE4_headscans_dur,10);
        title(['apoE4 S' num2str(s)])
        xlabel('Headscan Duration [s]')
        ylabel('Count')
    end
end