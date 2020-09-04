function [spikes, sampleCheck] = UMS2000(varargin)

    defaults = {...
        'direction', 'up';...
        'probe', 'P2';...
        'trodeIndex', 1;...
        'trode', 1:16;...
        'Fs', 32000;...
        'filepath', '';...
        };
    
    [s, ~] = parse_args(defaults, varargin{:});
    
    if isempty(s.filepath)
        s.filepath = uigetdir();
        if isempty(s.filepath)
            spikes = [];
            sampleCheck = [];
            return
        end
    end
    s.filepath = fullfile([s.filepath filesep]);    
    switch s.direction
        case 'up'
            inv = -1;
        case 'down'
            inv = 1;
    end
    
    switch s.probe
        case 'P2'
            % for P-2 probe using neurolynx adapter
            channelOrder = [16 5 12 4 10 2 9 1 7 3 8 6 14 13 15 11 17 28 21 29 23 31 24 32 26 30 25 27 19 20 18 22];
        case 'E2'           
            % for E-2 probe:
            channelOrder = [16 12 10 9 7 5 4 2 15 11 14 13 8 6 3 1 17 21 23 24 26 28 29 31 18 22 19 20 25 27 30 32];
        otherwise
            channelOrder = 1:32;
    end
    
    s.channelOrder = channelOrder;    
    
% testing probe site screening for light effect....
nlxcsc2mat2(s.filepath,'Channels','Events');
load(fullfile(s.filepath, 'Events.mat'));


% 
% % event timestamps are in seconds, CSC timestamps are in microseconds
startRecording = find(strcmp(Events_EventStrings, 'Starting Recording'), 1, 'first');
startRecording = Events_TimeStamps(startRecording);
% stopRecording = find(strcmp(Events_EventStrings, 'Stopping Recording'), 1, 'last');
% stopRecordingTime = Events_TimeStamps(startRecording);
% how many microseconds after recording start?



[Timestamps_all_pre, Fs_nlx, nValidSamples] = Nlx2MatCSC(fullfile(s.filepath, 'CSC1.ncs'),[1 0 1 1 0],0,1);% get # total sample pieces (output is nSamplePieces x 512)
if s.Fs ~= Fs_nlx
    error('wrong sample rate');
end
Timestamps_all_pre = Timestamps_all_pre * 1e-6;
Timestamps_all = repmat(Timestamps_all_pre, 512, 1);
Timestamps_all = Timestamps_all + ((0:511)/s.Fs)';
sampleCheck.nlx = sum(nValidSamples);
nValidSamples = numel(nValidSamples);
chunkSize = 1e5; % load spikes piecemeal so as not to overrun memory
nChunks = ceil(nValidSamples / chunkSize);
sampleRange = ((0:nChunks) * chunkSize);


trodeSize = length(s.trode); % formerly 4


sampleCheck.perChunk = zeros(nChunks, 1);
sampleCheck.range = zeros(nChunks, 2);
spikes = ss_default_params(s.Fs, {'display', 'trial_spacing'}, 0); % don't pad time between chunks/trials, see ss_detect
s.sampleCheck = sampleCheck;
spikes.nlx_settings = s; % special field 

h = waitbar(0, 'Detecting Spikes');  
for counter = 1:nChunks
    theseSamples = [(sampleRange(counter) + 1) min(sampleRange(counter + 1), nValidSamples)];
    sampleCheck.range(counter, :) = theseSamples;
    for ttcounter = 1:trodeSize
        cfilename = sprintf('CSC%d.ncs', channelOrder(s.trode(ttcounter)));
        [Timestamps, Samples, header] = Nlx2MatCSC(fullfile(s.filepath, cfilename),[1 0 0 0 1],1,2, theseSamples); %, [startRecording TrialStart_nlx(1)]);        
        display(['loaded ' cfilename]);
        % get length of data chunk (last chunk is shorter), bitVolts
        % conversion factor, initialize data array for spike detection
        Samples = reshape(Samples, numel(Samples), 1);        
        if ttcounter == 1
            nSamples = length(Samples);
            ADBitVolts = sscanf(header{16}, '-ADBitVolts %f');            
            data = zeros(1, nSamples, 4);
        end
        sampleCheck.perChunk(counter) = numel(Samples);

        
        Samples = Samples * ADBitVolts * 1e6 * inv;
        

        passBand = [300 8000]; % Hz
        [b, a] = butter(5, passBand/(s.Fs/2));
        Samples = filtfilt(b,a,Samples);
        data(1, :,ttcounter) = Samples;
        waitbar(((counter - 1) * trodeSize + ttcounter) / (nChunks * trodeSize));

    end   
    spikes = ss_detect(data, spikes);
    disp('spikes detected');
end
close(h);
% cfilename = sprintf('CSC%d.ncs', 5);
% [Samples, header] = Nlx2MatCSC(fullfile(filepath, cfilename),[0 0 0 0 1],1,1);
% Samples = reshape(Samples, numel(Samples), 1);
% sampleCheck.samples = Samples;

sampleCheck.chunksSummed = sum(sampleCheck.perChunk);
spikes.startRecording = startRecording;
% get sample indices of spikes, use them query nlx spiketimes
spikeIndices = round(spikes.unwrapped_times * spikes.params.Fs);
spikes.nlx_times = Timestamps_all(spikeIndices);

spikes = ss_align(spikes);
spikes = ss_kmeans(spikes);
spikes = ss_energy(spikes);
spikes = ss_aggregate(spikes);

splitmerge_tool(spikes, 'all', [], ['TT' num2str(s.trodeIndex)])
