function [hr, rr] = hr_rr(clrAvg,vidRate)
%#codegen
% Panic Button Heart and Respiratory Rate Function

% Inputs: 
% vidFrames - height x width x 3 x numFrames array, where each height x
%             width x 3 element is one video frame
% vidRate   - frame rate of the video
% numFrames - number of frames contained in the video
% Height    - height of the video frame (pixels)
% Width     - width of the video frame (pixels)

% Outputs:
% hr        - heart rate in Hz (beats per second)
% rr        - respiratory rate in Hz (breaths per second)

% Functional version of hr_rr_estimation.m that enables conversion to c.
% Written by Ryan S. McGinnis - ryan.mcginis14@gmail.com - June 15, 2016

% Copyright (C) 2016  Ryan S. McGinnis


% Calculate Heart Rate
% Apply some simple filtering
%TODO: update to custom filter definition so that we can use adaptive filter: cut_off=[0.5, 12]/(vidRate/2)
cut_off = [0.0334, 0.8005]; 
n=4;
[b,a] = butter(n,cut_off,'bandpass');
clrAvgFlt = filter(b,a,clrAvg);


% Take PCA of average color intensity time series
[~,score,~] = pca(clrAvgFlt,'NumComponents',1);


% Smooth iPPG waveform 
hr_wave = movingmean_v2(score,15,1);


% Estimate power spectrum with Welch's method
[pxx, freqs] = pwelch_v2(hr_wave,vidRate); 


% Calculate heart rate based on dominant frequency component
[~,max_ind] = max(pxx);
hr = freqs(max_ind);


% Calculate Respiratory Rate
% Calculate power spectrum of iPPG signal in sliding window to get an indication of HRV 
window = round(2 * vidRate);
overlap = round(window-1);
nfft = 2048;

[S,F,T] = sliding_fft(hr_wave,window,overlap,nfft,vidRate);

% Use windowed power spectra to estimate instantaneous heart rate
HR = zeros(size(S,2),1);
for col_ind = 1:size(S,2)
    [~,max_ind] = max(abs(S(:,col_ind)));
    HR(col_ind) = F(max_ind);
end

% Filter HR timeseries 
%TODO: update to custom filter definition so that we can use adaptive filter: cut_off=[0.05, hr*0.75]/(fs_RR/2)
fs_RR = 1/mean(diff(T));
cut_off = [0.0033, 0.0334];
n=4;
[b,a] = butter(n,cut_off,'bandpass');

HR_filt = filter(b,a,HR);


% Estimate power spectrum with Welch's method
[pxx_rr, freqs_rr] = pwelch_v2(HR_filt,fs_RR); %pwelch with default parameter values


% Calculate respiratory rate based on dominant frequency component
[~,max_ind] = max(pxx_rr);
rr = freqs_rr(max_ind);

end