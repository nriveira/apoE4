function [pss, filt_pss, num_headscans] = adCircTrack_x_findHeadscans(coords)
% ADCIRCTRACK_X_FINDHEADSCANS
%   Given a set of coordinates from the LED tracking, characterize headscans as in Monaco et al.
    % First, filter the position data by median filtering it to remove jump
    % discontinuities
    TRACK_RADIUS = 45;

    t = coords(:,1);
    fps = round(1/median(diff(t)));

    x = movmedian(coords(:,2), floor(fps/2));
    y = movmedian(coords(:,3), floor(fps/2));
    
    % Get polar coordinate transformation
    polar_angle = (pi+atan2(y, x))*TRACK_RADIUS;
    polar_r = sqrt(x.^2 + y.^2) - TRACK_RADIUS;

    dadt = diff(polar_angle).*fps;
    drdt = diff(polar_r);

    % Filter out any jumps over 30 cm/s to account for new laps
    dadt(dadt > 50) = nan;
    dadt = movmedian(dadt, fps, 'omitnan');
       
    % Variables to be changed later
    RADIAL_THRESH = 2;
    RUNNING_THRESH = 6;

    % Knierem algorithm for headscanning (Using 4 second buffer)
    running_buffer = zeros(fps*4, 1);
    pss = zeros(length(x), 1);
    buffer_count = 0;

    % At each frame, check conditions (Making sure to convert thresh)
    for i = 1:length(dadt)
        % Fill running buffer first
        % Make sure the running thresh and dadt are same units!
        if(abs(dadt(i)) > RUNNING_THRESH)
            running_buffer(1:end-1) = running_buffer(2:end);
            running_buffer(end) = i;
            buffer_count = buffer_count+1;
    
        % If the buffer has values, then continue looking for headscans
        elseif(buffer_count >= length(running_buffer))
            radial_buffer = polar_r(running_buffer);
            drdt_buffer = drdt(running_buffer);
    
            % If the polar coordinates are above radial threshold (off the
            % track)
            pss_c1 = abs(polar_r(i)) > RADIAL_THRESH;
    
            % If the polar coordinates are out of IQR range
            pss_c2 = polar_r(i) < prctile(radial_buffer,25) & polar_r(i) >= prctile(radial_buffer,75);
    
            % If drdt is outside of IQR range
            pss_c3 = drdt(i) < prctile(drdt_buffer,25) & drdt(i) >= prctile(drdt_buffer,75);
            
            if(pss_c1 || pss_c2 || pss_c3)
                pss(i+1) = 1;
            end
        end
    end
    
%     figure(1); clf; hold on;
%     plot(x, y, '.b')
%     plot(x(pss==1), y(pss==1), '.k');

    % Filter PSS using further criteria
    filt_pss = zeros(length(pss), 2);
    min_duration = fps*0.4; % In seconds
    min_size = 2.5; % In cm
    start_index = 0;
    prev = 0;
    filt_ind = 1;
    
    for i = 1:length(pss)
        if(pss(i) == 1)
            if(prev == 0)
                start_index = i;
            end
        else
            % Check conditions for min duration and min length
            if(prev == 1 && (i-start_index >= min_duration) && sumabs(drdt(start_index:i)) >= min_size)
                filt_pss(filt_ind, :) = [start_index, i];
                filt_ind = filt_ind+1;
            end
        end
        prev = pss(i);
    end
    filt_pss(filt_ind:end,:) = [];
    num_headscans = filt_ind;
end