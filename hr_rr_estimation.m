%% The Panic Button HR and RR estimation  
% This script uses a frequency-domain approach to estimate heart and
% respiratory rates from iPhone video of the finger tip.  More details
% about the method are provided here: https://experiment.com/u/BL7ghQ
% Written by Ryan S. McGinnis - ryan.mcginis14@gmail.com - June 15, 2016

% Copyright (C) 2016  Ryan S. McGinnis
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% Load video from specified file
filename = 'hr_vid.mp4'; % Add filename and path for video file here
vidObj = VideoReader(filename);
vidRate = vidObj.FrameRate;
numFrames = vidObj.NumberOfFrames;
vidFrames = read(vidObj);


%% Define region of interest in video frames
% TODO: Add Temporal Difference method from Rong-Chao, et al. 2016  
roi = [1, vidObj.Height, 1, vidObj.Width];


%% Extract average color intensities from each frame
clrAvg = zeros(numFrames,3); %[R, G, B]
for frameIndex=1:numFrames
    % Extract RGB video frame
    vidFrame = vidFrames(:,:,:,frameIndex);

    % Average colors across ROI and write to array
    for colorIndex=1:3
        clrAvg(frameIndex,colorIndex) = mean(mean(vidFrame(roi(1):roi(2),roi(3):roi(4),colorIndex)));
    end
end


% Clean up memory
clear vidFrames


%% Calculate Heart Rate
% Apply some simple filtering
cut_off = [0.5, 12];
n=4;
[b,a] = butter(n,cut_off/(vidRate/2),'bandpass');
clrAvgFlt = filtfilt(b,a,clrAvg);


% Take PCA of average color intensity time series
[~,score,~] = pca(clrAvgFlt,'NumComponents',1);


% Smooth iPPG waveform and plot
hr_wave = smooth(score,15,'moving');

figure;
set(gcf,'name','iPPG waveform raw (blue) and filtered (red)');
hold on;
plot(ts,score);
plot(ts,hr_wave,'r');
xlabel('Time (s)'); ylabel('HR Waveform');


% Estimate power spectrum with Welch's method and plot
window = [];
overlap = [];
nfft = [];
[pxx, freqs] = pwelch(hr_wave,window,overlap,nfft,vidRate); %pwelch with default parameter values

figure;
set(gcf,'name','iPPG power spectrum - peak corresponds to heart rate');
plot(freqs,10*log10(pxx));
ylabel('dB');
xlabel('Frequency (Hz)');
xlim([0, 220/60]); % set limit to physiological hr limits


% Calculate heart rate based on dominant frequency component
[~,max_ind] = max(pxx);
hr = freqs(max_ind);


%% Calculate Respiratory Rate
% Calculate spectrogram of iPPG signal to get an indication of HRV and plot
window = round(2 * vidRate);
overlap = round(window-1);
nfft = 2048;

figure;
set(gcf,'name','iPPG spectrogram');
spectrogram(hr_wave,window,overlap,nfft,vidRate,'yaxis');


% Calculate spectrogram again, but output power spectrum for each window 
[S,F,T] = spectrogram(hr_wave,window,overlap,nfft,vidRate,'yaxis');


% Use windowed power spectra to estimate instantaneous heart rate
HR = zeros(1,size(S,2));
for col_ind = 1:size(S,2)
    [~,max_ind] = max(abs(S(:,col_ind)));
    HR(col_ind) = F(max_ind);
end


% Filter HR timeseries and plot
fs_RR = 1/mean(diff(T));
cut_off = [0.05, hr*0.75];
n=4;
[b,a] = butter(n,cut_off/(fs_RR/2),'bandpass');

HR_filt = filtfilt(b,a,HR);

figure; 
set(gcf,'name','Instantaneous heart rate raw (blue) and filtered (red)'); 
hold on;
plot(T, HR * 60)
plot(T, HR_filt * 60, 'r');
ylabel('HR (bpm)');
xlabel('Time (s)');


% Estimate power spectrum with Welch's method and plot
window = [];
overlap = [];
nfft = [];
[pxx_rr, freqs_rr] = pwelch(HR_filt,window,overlap,nfft,fs_RR); %pwelch with default parameter values

figure;
set(gcf,'name','Heart rate power spectrum - peak corresponds to respiratory rate');
plot(freqs_rr,10*log10(pxx_rr));
ylabel('dB');
xlabel('Frequency (Hz)');
xlim([0, hr]);


% Calculate respiratory rate based on dominant frequency component
[~,max_ind] = max(pxx_rr);
rr = freqs_rr(max_ind);


%% Print results to the command line
fprintf('Data Length: %u seconds\n',round(numFrames / vidRate));
fprintf('HR: %f beats/min\n',hr * 60); 
fprintf('RR: %f breaths/min\n',rr * 60); 



