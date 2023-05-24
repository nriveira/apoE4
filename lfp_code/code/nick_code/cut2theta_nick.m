function [time, W, kInd] = cut2theta_nick(theta, Fs);

%  note from Laura:  input "theta" is the raw EEG trace.
%  output "W" is the cut theta cycles.
% Edited by Nick (5/23/2023) to output relative time to add spikes to wave.

%Fs = 1893;
Ftheta = 8;    % Hz
Ftol = 0.005;  % s

kstart = 1 - floor(0.3*Fs);
kstop =  1 + floor(0.3*Fs);
time = (1:kstop-kstart+1)/Fs - 0.3;

W = zeros(length(time),10000);

thetaBP = fftbandpass(theta,Fs,Ftheta-5,Ftheta-4,Ftheta+4,Ftheta+5);

kStart = [];
kStop = [];
kstart = 0;
j = 0;

for k=2:length(theta)-1
    % if thetaBP(k-1) > 0 & thetaBP(k) < 0 % downward slope
    %if thetaBP(k-1) < 0 & thetaBP(k) > 0 %upward slope
    if (thetaBP(k-1) > thetaBP(k))  & (thetaBP(k+1) > thetaBP(k)) %trough 
    %if (thetaBP(k-1) < thetaBP(k))  & (thetaBP(k+1) < thetaBP(k))  %peak
        kstartold = kstart;
        kstart = k - floor(0.3*Fs);
        kstop =  k + floor(0.3*Fs);
        time = (1:kstop-kstart+1)/Fs - 0.3;
        if kstart > 0 & kstop < length(thetaBP)
                dk = kstart - kstartold;
                kStart = [kStart; kstart];
                kStop = [kStop; kstop];
                if dk < Fs/Ftheta + Fs*Ftol & dk > Fs/Ftheta - Fs*Ftol
                    j = j + 1;
                    %Wd(:,j) = diff(theta(kstart:kstop));   
                    W(:,j) = theta(kstart:kstop);  
                end
        end
    end
    
end

W = W(:,1:j);
kInd = [kStart kStop];
