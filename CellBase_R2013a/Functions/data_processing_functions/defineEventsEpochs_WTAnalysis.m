function [events,epochs] = defineEventsEpochs_WTAnalysis
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
epochs(i,:) = {'WaitingTimeEndRate1',    'ResponseEnd',       [-1.5 -1],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'WaitingTimeEndRate2',    'ResponseEnd',       [-1 -0.3],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'WaitingTimeEndRate3',    'ResponseEnd',       [-0.3 0],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'WaitingTimeEndRate4',    'ResponseEnd',       [0 0.3],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'PostStart',    'ResponseStart',       [0 1],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart2',    'ResponseStart',       [0.75 1.25],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart3',    'ResponseStart',       [1.25 2],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart4',    'ResponseStart',       [2 3],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'WaitingTimeEndRate5',    'ResponseEnd',       [-2 -1.5],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'WaitingTimeEndRate6',    'ResponseEnd',       [-2.5 -2],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'WaitingTimeEndRate6',    'ResponseEnd',       [-3 -2.5],       'WaitingRewardTime'};    i = i + 1;
epochs(i,:) = {'PostStart5',    'ResponseStart',       [1.25 1.75],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart6',    'ResponseStart',       [1.75 2.25],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart7',    'ResponseStart',       [2.25 2.75],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart8',    'ResponseStart',       [2.75 3.25],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart9',    'ResponseStart',       [3.25 3.75],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart10',    'ResponseStart',       [3.75 4.25],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart11',    'ResponseStart',       [4.25 4.75],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart12',    'ResponseStart',       [4.75 5.25],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart13',    'ResponseStart',       [5.25 5.75],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart14',    'ResponseStart',       [5.75 6.25],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart15',    'ResponseStart',       [5.75 6.25],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart16',    'ResponseStart',       [6.25 6.75],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart17',    'ResponseStart',       [6.75 7.25],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart18',    'ResponseStart',       [7.25 7.75],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStart19',    'ResponseStart',       [7.75 8.25],     'MovementTime'};    i = i + 1;
epochs(i,:) = {'PostStartminus1',    'ResponseStart',       [-0.3 0.75],     'MovementTime'};    i = i + 1;





