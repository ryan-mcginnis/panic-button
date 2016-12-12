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
% hr        - heart rate in BPM (beats per minute)
% rr        - respiratory rate in BPM (breaths per minute)

% Functional version of hr_rr_estimation.m that enables conversion to c.
% Written by Ryan S. McGinnis - ryan.mcginis14@gmail.com - June 15, 2016

% Copyright (C) 2016  Ryan S. McGinnis


% Calculate Heart Rate
% Apply some simple filtering
cut_off = 2*[0.5, 12]/vidRate;
n=4;
[b,a] = butter_bp(n,cut_off);
clrAvgFlt = filter(b,a,clrAvg - (clrAvg(:,1).^0)*clrAvg(1,:));


% Take PCA of average color intensity time series
[~,score,~] = pca(clrAvgFlt,'NumComponents',1);


% Smooth iPPG waveform 
hr_wave = movingmean_v2(score,15,1);


% Estimate power spectrum with Welch's method
[pxx, freqs] = pwelch_v2(hr_wave,vidRate); 


% Find peaks and calculate instantaneous hr
hr_wave = real(hr_wave);
[pks,locs] = findpeaks(hr_wave,'MinPeakDistance',0.27*vidRate,'MinPeakHeight',0); %find peaks that are seperated by at least 1/(220 bpm / 60 s/m)
hrs = 1 ./ (diff(locs) / vidRate); % bps
hrs = hrs( hrs <= 220/60 & hrs >= 15/60);


%Estimate HR with recursive bayes
pxx = recursive_bayes(freqs, pxx, hrs);


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
fs_RR = 1/mean(diff(T));
if hr~=0
    cut_off = 2*[0.05, hr*0.75]/fs_RR;
else
    cut_off = 2*[0.05, 0.75]/fs_RR; %assume hr is actualy 60 bpm
end
n=4;
[b,a] = butter_bp(n,cut_off);

HR_filt = filter(b,a,HR-HR(1));


% Estimate power spectrum with Welch's method
[pxx_rr, freqs_rr] = pwelch_v2(HR_filt,fs_RR); %pwelch with default parameter values


% Calculate respiratory rate based on dominant frequency component
[~,max_ind] = max(pxx_rr);
rr = freqs_rr(max_ind);


% Concert hr and rr to BPM (beats/breaths per minute)
hr = 60 * hr;
rr = 60 * rr;

end