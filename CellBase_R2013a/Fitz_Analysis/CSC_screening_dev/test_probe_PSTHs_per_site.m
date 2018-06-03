% testing probe site screening for light effect....


[fname, pname] = uiputfile('path', 'Choose CSC path...');
filepath = pname;

load(fullfile(filepath, 'Events.mat')); 
TrialStart_nlx = getBehaviorStartTimes(Events_Nttls, Events_EventStrings, Events_TimeStamps);

% event timestamps are in seconds, CSC timestamps are in microseconds
startRecording = find(strcmp(Events_EventStrings, 'Starting Recording'), 1, 'last');
startRecordingTime = Events_TimeStamps(startRecording);
% how many microseconds after recording start?
timeStampRange = [startRecordingTime TrialStart_nlx] * 1e6; % convert to microseconds

trialDuration = 11;


nChannels = 16;
h = waitbar(0, 'Loading CSCs');   
data = struct(...
    'data', [],...
    'avg', [],...
    'std', [],...
    'thresh', []...
    );
for channel = 1:nChannels
    cfilename = sprintf('CSC%d.ncs', channel);
%     TimeStamps, ChannelNumbers, SampleFrequencies, NumberValidSamples, Samples and NlxHeader 
    [TimeStamps, Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[1 0 0 0 1],1,4, timeStampRange);%, [startRecording TrialStart_nlx(1)]);
    Samples = reshape(Samples, numel(Samples), 1);
    ADBitVolts = sscanf(header{16}, '-ADBitVolts %f');
    Samples = Samples * ADBitVolts * 1e6;
%     % downsample to 20KHz
%     Samples = resample(Samples, 32, 20);
    
    if channel == 1
        data.data = zeros(numel(Samples), nChannels);
        data2 = data.data;
    end
    passBand = [300 8000]; % Hz
    [b, a] = butter(5, passBand/32000 * 2);
    filtData = filtfilt(b,a,Samples);
    

    
    
    data.data(:, channel) = filtData;
    
    % compare to JRCLUST method of using Savitzky-Golay differentiation
    % filter
    [b, g] = sgolay(1,9);
    filtData = conv(Samples, -32000 * g(:,2), 'same');
    data2(:, channel) = filtData;
    waitbar(channel/nChannels);
end
close(h);


%% detect spikes, calculate SD per channel
data.avg = mean(data.data);
data.std = std(data.data);
data.thresh = data.std * -2;

data.spikeTimes = cell(nChannels, 1);
for channel = 1:nChannels
    crossings = find((data.data(1:end-1,channel) > data.thresh(channel) &  data.data(2:end,channel) <= data.thresh(channel)))';
%     crossings = find((data.data(1:end-1,channel) > data.thresh(channel) &  data.data(2:end,channel) <= data.thresh(channel)) |...
%         (data.data(1:end-1,channel) < -1 * data.thresh(channel) &  data.data(2:end,channel) >= -1 * data.thresh(channel)))';
    data.spikeTimes{channel} = (crossings / 32000) + startRecordingTime;
end

%% for each channel, calculate laser PSTH
PortID = eventPortFromEventStrings(Events_EventStrings);
Pulses = PortID == 0 & Events_Nttls == 128 & Events_TimeStamps < TrialStart_nlx(1);
PulseTimes = Events_TimeStamps(Pulses);

% 10Hz stimulation, 10ms bins
edges  = 0:.01:0.1;
histCountsByChannel = zeros(length(edges) - 1, nChannels);

for counter = 1:length(PulseTimes)
    pulseTime = PulseTimes(counter);
    for channel = 1:nChannels
        spikesZeroed = data.spikeTimes{channel} - pulseTime;
        histCountsByChannel(:,channel) = histCountsByChannel(:,channel) + histcounts(spikesZeroed, edges)';        
    end
end

histCountsByChannel_norm = (histCountsByChannel - mean(histCountsByChannel(5:end, :))) ./ mean(histCountsByChannel(5:end, :));
    




