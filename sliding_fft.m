function [S,F,T] = sliding_fft(data,window,overlap,nfft,fs)

S = zeros(round(nfft/2)+1, round((length(data)-window)/(window-overlap)));
data_len = length(data);

s = 1;
e = window;
i=1;
while e <= data_len
    % extract window of data
    temp = data(s:e,:);

    % calculate FFT 
    temp = fft(temp, nfft);

    % step 4: calculate the "periodogram" by taking the absolute value squared
    temp = abs(temp).^2;

    % save the results in the storage variable
    S(:, i) = temp(1:round(nfft/2)+1,1);

    % iterate indices
    i = i + 1;
    s = s + (window-overlap);
    e = s + window;
end

F = fs*(0:(nfft/2)).'/nfft;
T = (0:size(S,2)-1)/fs;

end