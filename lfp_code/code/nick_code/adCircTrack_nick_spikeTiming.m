function group = adCircTrack_nick_spikeTiming(group)
    Fs = 2000;
    tic
%ADCIRCTRACK_NICK_SPIKETRACK Summary of this function goes here
%   Detailed explanation goes here
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\lfp_e3e4';
    boxplot_ind = 1;

    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                tetNums = group(g).rat(r).day(d).tetNums;
                tt = group(g).rat(r).day(d).thetaTet;
                cscFn = ['CSC' num2str(tt) '.ncs'];

                phase_index = 0;
                thetaPhase_spkTms = zeros(1, 90);
                slowGammaPhase_spkTms = zeros(1, 90);
                fastGammaPhase_spkTms = zeros(1, 90);

                theta_slowGamma_spkTms = zeros(12, 12);
                theta_fastGamma_spkTms = zeros(12, 12);
                slowGamma_fastGamma_spkTms = zeros(12, 12);

                average_theta = [];
                average_slowGamma = [];
                average_fastGamma = [];

                for b = 1:length(group(g).rat(r).day(d).begin)
                    cd(group(g).rat(r).day(d).begin(b).dir)
                    % Load in lfp data and eeg indices
                    lfpStruct = read_in_lfp(cscFn);
                    
                    ts = lfpStruct.ts;
                    bpLFP = bandpass(zscore(lfpStruct.data), [0.5 200], lfpStruct.Fs);
                    unit = group(g).rat(r).day(d).begin(b).unit;

                    unitInd = reshape([unit.ID],2,[]);
                    unitInd(2,:) = [];
                    unit = unit(unitInd == tt);

                    if(b == 1) % Initialize once but get all begins
                        spikes_per_unit = [];
                        for u = 1:length(unit)
                            spikes_per_unit(u).ID = unit(u).ID(2);
                            spikes_per_unit(u).spikes = [];
                        end
                    end

                    for i = 1:length(group(g).rat(r).day(d).begin(b).eegInds)
                        eegStart = group(g).rat(r).day(d).begin(b).eegInds(i,1);
                        eegStop = group(g).rat(r).day(d).begin(b).eegInds(i,2);

                        startTime = ts(eegStart);
                        stopTime = ts(eegStop);          

                        % Run cut2theta to get the timings to align spikes
                        % to
                        eegSig = bpLFP(eegStart:eegStop);
                        thetaBP = fftbandpass(eegSig,Fs,3,4,12,13);
                        slowGammaBP = fftbandpass(eegSig, Fs, 20, 25, 45, 50);
                        fastGammaBP = fftbandpass(eegSig, Fs, 60, 65, 140, 145);

                        % Make the phase indices
                        theta_phase = rad2deg(angle(hilbert(thetaBP)))+180;
                        slowGamma_phase = rad2deg(angle(hilbert(slowGammaBP)))+180;
                        fastGamma_phase = rad2deg(angle(hilbert(fastGammaBP)))+180;

                        average_theta = [average_theta, Fs/(2*pi)*median(diff(unwrap(angle(hilbert(thetaBP)))))];
                        average_slowGamma = [average_slowGamma, Fs/(2*pi)*median(diff(unwrap(angle(hilbert(slowGammaBP)))))];
                        average_fastGamma = [average_fastGamma, Fs/(2*pi)*median(diff(unwrap(angle(hilbert(fastGammaBP)))))];
                       
                        %figure(1); clf; plot(eegSig, 'DisplayName','EEG Signal'); hold on; legend();

                        % Add spikes to the times given
                        for u = 1:length(unit)
                            spkTms = unit(u).spkTms([unit(u).spkTms] > startTime & [unit(u).spkTms] < stopTime);
                            spkTms = floor((spkTms - startTime)*lfpStruct.Fs)+1;

                            % Spike normalization factor for histogram
                            % (1/frequency)
                            spike_rate = 1/length(unit(u).spkTms);
                            
%                             plot(spkTms, eegSig(spkTms), '*', 'DisplayName',['Unit ' num2str(u)])

                            % Fill in spike times for each band
                            thetaPhase_spkTms = thetaPhase_spkTms + (spike_rate*histcounts(theta_phase(spkTms), 90, 'BinLimits', [0, 360]));
                            slowGammaPhase_spkTms = slowGammaPhase_spkTms + (spike_rate*histcounts(slowGamma_phase(spkTms), 90, 'BinLimits', [0, 360]));
                            fastGammaPhase_spkTms = fastGammaPhase_spkTms + (spike_rate*histcounts(fastGamma_phase(spkTms), 90, 'BinLimits', [0, 360]));
                    
                            % Get the indices for the matrix
                            theta_ind = floor(theta_phase(spkTms)/30)+1;
                            sg_ind = floor(slowGamma_phase(spkTms)/30)+1;
                            fg_ind = floor(fastGamma_phase(spkTms)/30)+1;
                    
                            % Fill in the matrix using the indices above
                            theta_slowGamma_spkTms(theta_ind, sg_ind) = theta_slowGamma_spkTms(theta_ind, sg_ind)+spike_rate;
                            theta_fastGamma_spkTms(theta_ind, fg_ind) = theta_fastGamma_spkTms(theta_ind, fg_ind)+spike_rate;
                            slowGamma_fastGamma_spkTms(sg_ind, fg_ind) = slowGamma_fastGamma_spkTms(sg_ind, fg_ind)+spike_rate;

                            % Also add spikes per unit info
                            spikes_per_unit(u).spikes = [spikes_per_unit(u).spikes, theta_phase(spkTms)];
                        end
                    end
                    boxplot_data(boxplot_ind).theta = average_theta;
                    boxplot_data(boxplot_ind).slowGamma = average_slowGamma;
                    boxplot_data(boxplot_ind).fastGamma = average_fastGamma;
                    boxplot_ind = boxplot_ind+1;
                end

                
                figure(g); hold on;
                subplot(3,4,d); 
                histogram('BinEdges', 0:4:720, 'BinCounts', [thetaPhase_spkTms, thetaPhase_spkTms]);
                if(d == 1)
                    title('Theta Phase');
                    xlabel('Phase (degrees)'); 
                    ylabel('Normalized Spike Rate'); 
                else
                    title(['Day ' num2str(d)])
                end

                subplot(3,4,d+4); 
                histogram('BinEdges', 0:4:720, 'BinCounts', [slowGammaPhase_spkTms, slowGammaPhase_spkTms]); 
                if(d == 1)
                    title('Slow Gamma Phase');                    
                    xlabel('Phase (degrees)'); 
                    ylabel('Normalized Spike Rate'); 
                else
                    title(['Day ' num2str(d)])
                end

                subplot(3,4,d+8); 
                histogram('BinEdges', 0:4:720, 'BinCounts', [fastGammaPhase_spkTms, fastGammaPhase_spkTms]); 
                if(d == 1)
                    title('Fast Gamma Phase'); 
                    xlabel('Phase (degrees)'); 
                    ylabel('Normalized Spike Rate'); 
                else
                    title(['Day ' num2str(d)])
                end

                figure(g+2); hold on; 
                subplot(3,4,d);
                imagesc(1:30:360, 1:30:360, theta_slowGamma_spkTms); colorbar; 
                if(d == 1)
                    title([group(g).name ' Theta-Slow Gamma Phase']);
                    colormap hot; 
                    ylabel('Slow Gamma Phase'); 
                    xlabel('Theta Phase');
                else
                    title(['Day ' num2str(d)])
                end
                subplot(3,4,d+4); 
                imagesc(1:30:360, 1:30:360, theta_fastGamma_spkTms); colorbar; 
                if(d == 1)
                    title([group(g).name ' Theta-Fast Gamma Phase']);
                    colormap hot;
                    ylabel('Fast Gamma Phase'); 
                    xlabel('Theta');
                else
                    title(['Day ' num2str(d)])
                end

                subplot(3,4,d+8); 
                imagesc(1:30:360, 1:30:360, slowGamma_fastGamma_spkTms); colorbar; 
                if(d == 1)
                    title([group(g).name ' Slow Gamma-Fast Gamma Phase']);
                    colormap hot; 
                    ylabel('Fast Gamma Phase'); 
                    xlabel('Slow Gamma Phase');
                else
                    title(['Day ' num2str(d)])
                end

%                 figure(figure_index); figure_index=figure_index+1; clf;
%                 for u = 1:length(spikes_per_unit)
%                     subplot(length(spikes_per_unit), 2, 2*(u-1)+1)
%                     histogram(spikes_per_unit(u).spikes, 90); 
%                     title([group(g).name ' Day ' num2str(d) ' Unit ' num2str(u)]);
%                     subplot(length(spikes_per_unit), 2, 2*u);
%                     histogram(mod((diff(spikes_per_unit(u).spikes)+360),360), 90); 
%                     title(['Number of spikes: ' num2str(length(spikes_per_unit(u).spikes))]);
%                 end
                %saveas(gcf, [saveDir filesep '20230615_SpikePhasePerCell_Day' num2str(d) group(g).name '.png'])
            end
        end
        figure(g); saveas(gcf, [saveDir filesep '20230615_SpikePhases_Day' num2str(d) group(g).name '.png'])
        figure(g+2); saveas(gcf, [saveDir filesep '20230615_SpikePhasePhase_Day' num2str(d) group(g).name '.png'])
    end
    toc
end
