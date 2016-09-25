%% Main script for illustrating call to hr_rr function
% Written by Ryan S. McGinnis - ryan.mcginis14@gmail.com - Sep. 24, 2016

%% Load video from specified file
filename = 'hr_vid.mp4'; % Add filename and path for video file here
vidObj = VideoReader(filename);
vidRate = vidObj.FrameRate;
numFrames = round(vidObj.Duration * vidRate);
vidFrames = read(vidObj);

%% Extract and aggregate color information from video frames
% Extract average color intensities from each video frame
% Define region of interest in video frames
% TODO: Add Temporal Difference method from Rong-Chao, et al. 2016  
roi = [1, vidObj.Height, 1, vidObj.Width];


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


%% Calculate heart and respiratory rates
tic;
[hr, rr] = hr_rr(clrAvg,vidRate);
toc