function [time, W, kInd] = cut2freq_nick(wave, freq, Fs);

%  note from Laura:  input "theta" is the raw EEG trace.
%  output "W" is the cut theta cycles.
%
% Adapted from cut2theta by Nick to look for other frequencies

%Fs = 1893;
%freq = 8;    % Hz
Ftol = 0.005;  % s
outerWindow = 5;
innerWindow = 4;

windowTime = 1/freq;
time = -windowTime:1/Fs:windowTime;

W = [];

waveBP = fftbandpass(wave,Fs,freq-outerWindow,freq-innerWindow,freq+innerWindow,freq+outerWindow);

kStart = [];
kStop = [];
kstart = 0;
j = 0;

for k=2:length(wave)-1
    % if thetaBP(k-1) > 0 & thetaBP(k) < 0 % downward slope
    %if thetaBP(k-1) < 0 & thetaBP(k) > 0 %upward slope
    if (waveBP(k-1) > waveBP(k))  & (waveBP(k+1) > waveBP(k)) %trough 
    %if (thetaBP(k-1) < thetaBP(k))  & (thetaBP(k+1) < thetaBP(k))  %peak
        kstartold = kstart;
        kstart = k - floor(Fs/freq);
        kstop =  k + floor(Fs/freq);
        if kstart > 0 & kstop < length(waveBP)
                dk = kstart - kstartold;
                if dk < Fs/freq + Fs*Ftol & dk > Fs/freq - Fs*Ftol
                    kStart = [kStart; kstart];
                    kStop = [kStop; kstop];
                    j = j + 1;
                    %Wd(:,j) = diff(theta(kstart:kstop));   
                    W(:,j) = wave(kstart:kstop);  
                end
        end
    end
    
end

kInd = [kStart kStop];