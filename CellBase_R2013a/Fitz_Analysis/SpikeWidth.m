function [p2t] = SpikeWidth(cellid, varargin)

%% optional parameters, first set defaults
defaults = {...
    'Fs', 32000;...
    };

[s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings



sw = loadcb(cellid,'Waveforms');   % load spike waveforms



sw_mean = squeeze(nanmean(sw,1));
% find max channel
[~, chMax ] = max(range(sw_mean, 2));
% find peak (spikes are upside down)
[~, ixPeak] = min(sw_mean(chMax, :));
% find trough (spikes are upside down)
[~, ixTrough] = max(sw_mean(chMax, ixPeak:end));
ixTrough = ixTrough + ixPeak - 1;
p2t = NaN;
if ixTrough > ixPeak && ixTrough ~= size(sw_mean, 3)
    p2t = (ixTrough - ixPeak) * 1/s.Fs;
else
    p2t = NaN; % trough has to occur later than peak and earlier than end of waveform
end    
