function headscan_struct = analyze_headscans(fileLoc, parameters)
% ANALYZE HEADSCANS 
%   Run headscan analsis on output of DLC, called by plot_headscans
    headscan_struct.dlc_struct = get_dlc_struct(fileLoc);

    % Convert radial distance to cm
    radial = median([headscan_struct.dlc_struct(:).predicted_r],2);
    radial_cm = radial./(median([headscan_struct.dlc_struct(:).r_length])/parameters.track_in_cm);
    headscan_struct.radial_cm = abs(radial_cm - median(radial_cm, 'omitnan'));

    % Also convert radial distance velocities to cm/s
    drdt = diff(headscan_struct.radial_cm).*parameters.fps;
    headscan_struct.drdt = movmedian(drdt, parameters.median_window);
    
    % Get polar angle, then multiply by radius to get angular velocity
    polar_angle = median([headscan_struct.dlc_struct(:).polar_coord],2);
    dadt = (parameters.track_in_cm* parameters.fps) * median([headscan_struct.dlc_struct(:).dadt],2);
    dadt(dadt > parameters.v_max | dadt < -parameters.v_max) = nan;

    headscan_struct.angle = polar_angle;
    headscan_struct.dadt = movmedian(dadt, parameters.median_window, 'omitnan');
    headscan_struct.total_time = length(headscan_struct.dadt)/parameters.fps;
    
    % Amount of pause time in seconds
    headscan_struct.pause_time = sum(parameters.track_in_cm*headscan_struct.dadt > parameters.pause_thresh) / parameters.fps;

    % Count the number of laps using 32 evenly spaced intervals to track
    % the position
    track_spacing = 2*pi/32;
    current_angle = median(headscan_struct.angle(1:10), 'omitnan');
    %headscan_struct.laps = zeros(size(headscan_struct.angle));
    headscan_struct.laps = unwrap(median([headscan_struct.dlc_struct(:).polar_coord],2)) / (2*pi);

%     for pos = 1:length(headscan_struct.angle)
%        if((headscan_struct.angle(pos) >= mod(current_angle, 2*pi)) & (headscan_struct.angle(pos) <= mod(current_angle, 2*pi)+track_spacing))
%             current_angle = current_angle+track_spacing;
%        end
%        headscan_struct.laps(pos) = current_angle/(2*pi);
%     end

    % Knierem algorithm for headscanning (Using 4 second buffer)
    running_buffer = zeros(parameters.fps*4, 1);
    headscan_struct.pss = zeros(length(headscan_struct.dadt), 1);
    buffer_count = 0;

    % At each frame, check conditions (Making sure to convert thresh)
    for i = 1:length(headscan_struct.dadt)
        % Fill running buffer first
        % Make sure the running thresh and dadt are same units!
        if(parameters.track_in_cm*headscan_struct.dadt(i) > parameters.running_thresh)
            running_buffer(1:end-1) = running_buffer(2:end);
            running_buffer(end) = i;
            buffer_count = buffer_count+1;
    
        % If the buffer has values, then continue looking for headscans
        elseif(buffer_count >= length(running_buffer))
            radial_buffer = headscan_struct.radial_cm(running_buffer);
            drdt_buffer = headscan_struct.drdt(running_buffer);
    
            % If the polar coordinates are above radial threshold (off the
            % track)
            pss_c1 = headscan_struct.radial_cm(i) > parameters.radial_thresh/parameters.track_in_cm;
    
            % If the polar coordinates are out of IQR range
            pss_c2 = headscan_struct.radial_cm(i) < prctile(radial_buffer,25) & headscan_struct.radial_cm(i) >= prctile(radial_buffer,75);
    
            % If drdt is outside of IQR range
            pss_c3 = headscan_struct.drdt(i) < prctile(drdt_buffer,25) & headscan_struct.drdt(i) >= prctile(drdt_buffer,75);
            
            if(pss_c1 || pss_c2 || pss_c3)
                headscan_struct.pss(i) = 1;
            end
        end
    end

    % If there are events that are at least 0.4 seconds between them,
    % connect them
    duration_threshold = ones(floor(parameters.fps*0.4),1); % In seconds
    threshold_pss = conv(headscan_struct.pss, duration_threshold, "same");
    headscan_struct.pss = [threshold_pss >= 5];


    % Filter PSS using further criteria
    headscan_struct.filt_pss = zeros(length(headscan_struct.pss), 2);
    min_duration = parameters.fps*0.4; % In seconds
    min_size = 2.5; % In cm
    start_index = 0;
    prev = 0;
    filt_ind = 1;
    
    for i = 1:length(headscan_struct.pss)
        if(headscan_struct.pss(i) == 1)
            if(prev == 0)
                start_index = i;
            end
        else
            % Check conditions for min duration and min length
            if(prev == 1 && (i-start_index >= min_duration) && sumabs(headscan_struct.drdt(start_index:i)) >= min_size)
                headscan_struct.filt_pss(filt_ind, :) = [start_index, i];
                filt_ind = filt_ind+1;
            end
        end
        prev = headscan_struct.pss(i);
    end
    headscan_struct.filt_pss(filt_ind:end,:) = [];
    headscan_struct.num_headscans = filt_ind;
    headscan_struct.num_laps = range(headscan_struct.laps);
end