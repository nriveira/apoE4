% Plot individual head scan file
function output = plot_headscans(fileLoc, param)    
    headscan_struct = analyze_headscans(fileLoc, param);
    
    %% Plot Headscan variables
    t = (0:length(headscan_struct.dadt))./param.fps;
    figure(1); clf;
    subplot(3,2,1);
    plot(t, headscan_struct.radial_cm)
    title('Radial Distance from center')
    xlabel('Time [s]')
    ylabel('Position [cm]')
    
    subplot(3,2,2); 
    histogram(headscan_struct.radial_cm)
    
    subplot(3,2,3); hold on;
    plot(t(1:end-1), (param.fps*param.track_in_cm)*headscan_struct.dadt, 'b')
    title('Angular Velocity')
    xlabel('Time [s]')
    ylabel('Velocity [cm/s]')
    
    subplot(3,2,4); 
    histogram(headscan_struct.dadt)
    
    subplot(3,2,5);
    plot(t(1:end-1), headscan_struct.drdt)
    title('DfC Velocity')
    xlabel('Time [s]')
    ylabel('Velocity [cm/s]')
    
    subplot(3,2,6);
    histogram(headscan_struct.drdt)
    
    %% Plot the track and distance traveled
    figure(2); clf;
    subplot(2,2,1); hold on;
    plot(median([headscan_struct.dlc_struct(:).x_centered],2), median([headscan_struct.dlc_struct(:).y_centered],2), '.');
    title('Median points')
    axis square;
    
    subplot(2,2,2)
    x = median([headscan_struct.dlc_struct(:).x_centered],2);
    y = median([headscan_struct.dlc_struct(:).y_centered],2);

    plot(x(headscan_struct.pss==1), y(headscan_struct.pss==1), '.');
    title('Putative Scanning Sample')
    axis square;
    
    subplot(2,2,3);
    plot(t, headscan_struct.laps)
    title('Laps completed')

    subplot(2,2,4);
    plot(headscan_struct.pause_time, '*k')
    title('Pause Time')
    
    %% Plot individual events
    figure(3); clf; hold on;
    for i = 1:length(headscan_struct.filt_pss)
        x = median([headscan_struct.dlc_struct(:).x_centered],2); 
        y = median([headscan_struct.dlc_struct(:).y_centered],2);
        plot(x(headscan_struct.filt_pss(i,1):headscan_struct.filt_pss(i,2)), y(headscan_struct.filt_pss(i,1):headscan_struct.filt_pss(i,2)),'.');
    end
    viscircles([0,0], median([headscan_struct.dlc_struct(:).r_length]),'Color','k')
    title(['n=' num2str(length(headscan_struct.filt_pss))]);
    axis square
end