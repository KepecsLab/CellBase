function [p, I] = calcSALT(cellid, varargin)

%% optional parameters, first set defaults
defaults = {...
    'baselineEvent', 'BurstOn';... 
    'baselineWindow', [-0.5 0];...
    'testEvent', 'PulseOn';...
    'testWindow', [0 0.5];...
    'trials', [];...  % if empty, use all trials from StimEvents 
    'dt', 0.001; % resolution of bin raster in seconds
    'wn', 0.001; % window size for SALT calculation (see salt.m)
    };

[s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings



ST = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes for stimulation events
SE = loadcb(cellid,'StimEvents');


% get spike times for baseline and test events
pos_bl = findcellstr(ST.events(:,1),s.baselineEvent); % position of basleline event
pos_test = findcellstr(ST.events(:,1),s.testEvent);
if pos_bl == 0 || pos_test == 0
    error('Events not found for baseline and test periods');
end


if isempty(s.trials)
    s.trials = 1:length(SE.(s.baselineEvent));
end

blTrials = intersect(find(~isnan(SE.(s.baselineEvent))), s.trials);
testTrials = intersect(find(~isnan(SE.(s.testEvent))), s.trials);


% Calculate bin rasters
spiketimes_bl = stimes2binraster(ST.event_stimes{pos_bl}(blTrials),s.baselineWindow(1):s.dt:s.baselineWindow(2),s.dt);
spiketimes_test = stimes2binraster(ST.event_stimes{pos_test}(testTrials),s.testWindow(1):s.dt:s.testWindow(2),s.dt);

% SALT
[p I] = salt(spiketimes_bl,spiketimes_test,s.dt,s.wn);
