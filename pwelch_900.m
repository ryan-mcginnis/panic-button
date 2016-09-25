function [pxx, f] = pwelch_900(data, fs)
    %TODO update this to be hard coded so we don't calculate every time,
    %based on fixed data input size
    win_len = 512;
    num_windows = 12;    
    nfft = 2048;
    data_len = 900;
    overlap = 32; 

    % make a hamming window
    w = hamming(win_len);    
    
    % preallocate the space for the psd results
    mypsd = zeros(nfft, num_windows);
    
    % step 1: loop through the data, "win_len" points at a time, with "overlap" points overlap
    e=0;
    i=0;
    while e < data_len
        % step 2: apply a hamming window
        temp = data(1+overlap*i:win_len+i*overlap).*w;
       
        % step 3: calculate FFT and take just the first half
        temp = fft(temp, nfft);
        
        % step 4: calculate the "periodogram" by taking the absolute value squared
        temp = abs(temp).^2;
        
        % save the results in the storage variable
        mypsd(:, i+1) = temp;
        
        %update iterators
        e = win_len+(i+1)*overlap;
        i = i+1; 
    end
    
    % step 5: take the average of all the periodograms
    pxx = mean(mypsd,2);
    
    % throw away the 2nd half of mypsd
    pxx = pxx(1:round(nfft/2)+1);
    
    % normalizing factor
    pxx = pxx/(fs*sum(w.^2));
    
    % ignore the DC and Nyquist value
    pxx(2:end-1) = pxx(2:end-1) * 2;
    f = fs*(0:(nfft/2))/nfft;

end