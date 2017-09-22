function [events,epochs] = defineEventsEpochs_Dual2AFC
%DEFINEEVENTSEPOCHS_GONOGO   Define events and epochs for spike extraction.
%   [EVENTS,EPOCHS] = DEFINEEVENTSEPOCHS_GONOGO defines events and epochs
%   for spike extraction. 
%
%   EVENTS is a Nx4 cell array with columns corresponding to EventLabel,
%   EventTrigger1, EventTrigger2, Window. EventLabel is the name for
%   referencing the event. EventTrigger1 and EventTrigger2 are names of
%   TrialEvent variables (e.g. 'LeftPortIn'). For fixed windows, the two
%   events are the same; for variable windows, they correspond to the start
%   and end events. Window specifies time offsets relative to the events;
%   e.g. events(1,:) = {'OdorValveOn','OdorValveOn','OdorValveOn',[-3 3]};
%
%   EPOCH is a Nx4 cell array with columns corresponding to  EpochLabel, 
%   ReferenceEvent, Window, RealWindow. EventLabel is the name for 
%   referencing the epoch. ReferenceEvent should match an EventLabel in 
%   EVENTS (used for calculating the epoch rates). RealWindow is currently
%   not implemented (allocated for later versions).
%
%   DEFINEEVENTSEPOCHS_GONOGO defines events and epochs for auditory
%   go-nogo task.
%
%   See also MAKETRIALEVENTS2_GONOGO and DEFINEEVENTSEPOCHS_DEFAULT.

%   Edit log: BH 7/6/12 PM 7/03/14

% Define events and epochs
%              EventLabel       EventTrigger1      EventTrigger2      Window
i = 1;
events(i,:) = {'StimulusOnset',    'StimulusOnset',      'StimulusOnset',      [-6 6]};    i = i + 1;
events(i,:) = {'StimulusOffset',   'StimulusOffset',     'StimulusOffset',     [-6 6]};    i = i + 1;
events(i,:) = {'ResponseStart',    'ResponseStart',      'ResponseStart',      [-6 10]};    i = i + 1;
events(i,:) = {'ResponseEnd',      'ResponseEnd',        'ResponseEnd',          [-10 6]};    i = i + 1;




% Variable events
events(i,:) = {'StimulusSamplingDuration',    'StimulusOnset',        'StimulusOffset',     [-6 6]};    i = i + 1;
events(i,:) = {'WaitingRewardTime',                 'ResponseStart',        'ResponseEnd',     [-6 6]};    i = i + 1;
events(i,:) = {'MovementTime',                'StimulusOffset',       'ResponseStart',     [-6 6]};    i = i + 1;


% Define epochs for rate calculations
%               EpochLabel      ReferenceEvent      FixedWindow       RealWindow
i = 1;
epochs(i,:) = {'StimulusEndRate1',    'StimulusOffset',       [-0.3 -0.2],     'StimulusSamplingDuration'};    i = i + 1;
epochs(i,:) = {'StimulusEndRate2',    'StimulusOffset',       [-0.2 0.0],     'StimulusSamplingDuration'};    i = i + 1;
epochs(i,:) = {'StimulusEndRate3',    'StimulusOffset',       [0.1 0.3],     'StimulusSamplingDuration'};    i = i + 1;
epochs(i,:) = {'WaitingTimeEndRate3',    'ResponseEnd',       [-0.5 0],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'WaitingTimeEndRate4',    'ResponseEnd',       [0 0.5],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'PostStart',    'ResponseStart',       [0 0.75],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStartminus1',    'ResponseStart',       [-0.2 0.5],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart',    'ResponseStart',       [0 1.2],     'MovementTime'};    i = i + 1;





