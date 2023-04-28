function [onset,offset, ind] = DetectScanningEvent(post,posx,posy,circle)
% Detect scanning events on circular track. 
% The detection method followed Monaco et al., 2014 Nature Neuroscience

[scanning_ind,angle,r] = putativeSample(post,posx,posy,circle);
[onset, offset, ind] = ScanningEvent(scanning_ind,angle,r,post);

% figure; plot(posx,posy);hold on;
% for rr = 1: length(onset)
%     range = onset(rr):offset(rr);
%     plot(posx(range),posy(range),'r*')
% end

function [scanning_ind,angle,r] = putativeSample(post,posx,posy,circle)
r_thr = 5; %cm
speed_thr = 10/360*2*pi; %in radian 
buffer_length = 4; %sec

radius = circle.radius;
angle = atan2(posy,posx);
[speed_cir] = speedestimated(post,posx,posy,'circular',radius);

r = sqrt(posx.^2+posy.^2)-radius;
r_bar = diff(r)./diff(post);
r_bar(end+1)=0;

running_ind = speed_cir'>=speed_thr*radius & abs(r) < r_thr;
sampleL = round(1/mean(diff(post))*buffer_length);
if mod(sampleL,2) == 1; sampleL = sampleL-1; end
scanning_ind = false(size(posx));
for rr = sampleL/2:length(running_ind)-sampleL/2
    range = rr-sampleL/2+1:rr+sampleL/2;
    
    lr = prctile(r(range),25);lr_bar = prctile(r_bar(range),25);
    ur = prctile(r(range),75);ur_bar = prctile(r_bar(range),75);
    delta_r = ur-lr; delta_r_bar = ur_bar-lr_bar;
    
    %The three ways to be qualitfied as putative scanning samples
    temp1 = abs(r(range)) > r_thr;
    temp2 = r(range) < (lr-delta_r/2) | r(range) > (ur+delta_r/2);
    temp3 = r_bar(range) < (lr_bar-delta_r_bar/2) | r_bar(range) > (ur_bar+delta_r_bar/2);
    
    scanning_ind(range) = temp1 | temp2 | temp3;
end
scanning_ind(running_ind) = false;

function [onset, offset, ind] = ScanningEvent(scanning_ind,angle,r,post)
isi_thr = .4; % in sec
dur_thr = [.4 30]; % in sec
angle_thr = 45/360*2*pi; % maximum angle displacement
r_thr = 2.5; % in cm; radial displacement

sampfreq = 1/mean(diff(post));
onset = find(diff(scanning_ind)==1)+1;
offset = find(diff(scanning_ind)==-1);
%begin condition
if scanning_ind(1) == 1
    onset = [1 onset];        
end       
%end condition
if scanning_ind(end) == 1
    offset = [offset size(scanning_ind,2)];        
end

isi = onset(2:end)-offset(1:end-1);
isi_combine = isi < isi_thr*sampfreq;
onset([false isi_combine]) = [];
offset([isi_combine false]) = [];

dur = offset-onset;
onset(dur<dur_thr(1)*sampfreq | dur>dur_thr(2)*sampfreq) = [];
offset(dur<dur_thr(1)*sampfreq | dur>dur_thr(2)*sampfreq) = [];

rmv = false(size(onset));
for ii = 1:length(onset)
    range = onset(ii):offset(ii);
    angle_unwrap = unwrap(angle(range));
    if max(angle_unwrap)-min(angle_unwrap) > angle_thr
        rmv(ii) = true;
    end
    if max(r(range))-min(r(range)) < r_thr
        rmv(ii) = true;
    end
end
onset(rmv)=[];
offset(rmv)=[];

ind = false(size(scanning_ind));
for ii = 1:length(onset)
    range = onset(ii):offset(ii);
    ind(range) = true;
end
