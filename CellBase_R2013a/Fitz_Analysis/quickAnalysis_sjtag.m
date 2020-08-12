function quickAnalysis_sjtag(animalID,sessionID,sessionspec,protocoltag)
%QUICKANALYSIS2   Analysis of tetrode data.
%   QUICKANALYSIS2 is designed as an offline analysis tool for tetrode data
%   and behavior, which can be executed in an unsupervised manner on a
%   daily bases. It gives a quick overview of the experiment including
%   response profiles of clustered neurons, light-triggered PSTH and
%   psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   QUICKANALYSIS2(ANIMALNO,SESSIONID,SESSIONSPEC) performs the analysis
%   for a session specified by the first two input arguments. SESSIONSPEC
%   should be a 1x3 logical array indicating the need for behavior,
%   recording and stimulation analysis.
%ju
%   QUICKANALYSIS2(ANIMALNO,SESSIONID,SESSIONSPEC,PROTOCOLTAG) accepts a
%   PROTOCOLTAG argument to allow calls to trial event conversion programs
%   for different versions of the behavioral protocol.
%
%   See also ADDNEWCELLS, PREALIGNSPIKES and VIEWCELL2B.

% Input argument check
error(nargchk(0,4,nargin))

% Behavior, recording or both
if nargin < 4
    protocoltag = '';
end
if nargin < 3
    isbeh = 0;
    isrec = 1;
    isstim = 1;
else
    isbeh = sessionspec(1);
    isrec = sessionspec(2);
    isstim = sessionspec(3);
end

% Animal, session
% if nargin < 2
%     sessionID = '170909a';
% end
% if nargin < 1
%     animalID2 = 'nb046'; % Note FS- why two of these?
%     animalID = 'CD4';
% else
%     animalID2 = ['nb0' num2str(animalNO)];
%     animalID = ['n0' num2str(animalNO)];
% end

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

% Stop if error
dbstop if error



% Convert events file
if isrec
    nlxcsc2mat2(fullpth,'Channels','Events')
    if isbeh
%         TE = makeTE_CuedOutcome_Odor_Complete_Nlx(fullpth);
    end
end

% Create trial events structure
if isbeh
%     if isempty(protocoltag)
%     validTrials = filterTE(TE, 'reject', 0);
    
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter));%, 'reject', 0);
    end
    

    

end

% Update CellBase
if isrec
%     addnewcells('dir',[animalID filesep sessionID]) % change this back
    cellids = findcell('rat',animalID,'session',sessionID);
%     cellids = {'sj201_200808a_6.1'};
    disp(cellids)
end

% Response profiles
if isbeh && isrec
    % Prealign spikes for trial events
    problem_behav_cellid = [];
    for iC = 1:length(cellids)
        cellid = cellids(iC);
%         try
%             prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_CuedOutcome,'filetype','event','ifsave',1,'ifappend',0, 'writing_behavior', 'overwrite')
%         catch
%             disp('Error in prealignSpikes.');
%             problem_behav_cellid = [problem_behav_cellid cellid];
%         end
    end
    

end

% Light effects
if isrec && isstim

    % Create stimulus events
    MakeStimEvents_Bpod(fullpth,'PulseNttl',128, 'PulsePort', 0); % FS
    
    % Prealign spikes to stimulus events
    problem_stim_cellid = [];
    for iC = 1:length(cellids)
        cellid = cellids(iC);
        try
            prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim_sj,'filetype','stim','ifsave',1,'ifappend',0)
        catch
            disp('Error in prealignSpikes.');
            problem_stim_cellid = [problem_stim_cellid cellid];
        end
    end
    
    % View light-triggered raster and PSTH
    TrigEvent = 'Pulse';
    SEvent = 'Pulse';
    win = [-0.5 0.5];
    parts = 'all';
%     parts = '#BurstNPulse';
    dt = 0.001;
    sigma = 0.001;
    PSTHstd = 'on';
    ShEvent = {'Pulse'};
    ShEvColors = hsv(length(ShEvent{1}));
    ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
    for iCell = 1:length(cellids)
        continue;
        cellid = cellids(iCell);
        H = ensureFigure([cellids{iCell} '_laserStim'], 1);
        viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',ShEvColors,...
            'FigureNum',H,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
            'EventMarkerWidth',0,'PlotZeroLine','off')
%         maximize_figure(H)
        formatFigureCellbase;
        
        cellidt = cellid{1};
        cellidt(cellidt=='.') = '_';
        saveas(H, fullfile(fullpth, [cellidt '_LS.jpg']));
%         saveas(H, fullfile(fullpth, [cellidt '_LS.fig']));        
%         close(H)

        % research statement-  tagging figure
    end
end


% Cluster quality
% if isrec
    for iCell = 1:length(cellids)
        cellid = cellids{iCell};
        
        % get waveforms
        dummy1 = figure;
        dummy2 = figure;
        H = ensureFigure([cellid '_wave_acg'], 1);
        fighandle = [dummy1; dummy2; H];
%         st = loadcb(cellid, 'Spikes'); % spiketimes
        subplot(
        plotwaveforms(cellid, 'evoked', true, 'spont', true, 'compare', true, 'stim_period', [0 0.01], 'fighandle', repmat(H.Number, 3, 1));
        close(dummy1); close(dummy2);
        saveas(H, fullfile(fullpth, [cellid '_wave_acg.jpg']));
                       
    end
        
% end

% % Behavior
% if isbeh
%     auditory_gonogo_psychplot2(animalID,sessionID)
%     % auditory_gonogo_psychplot2(animalID,sessionID,[],(1:150))
%     H = gcf;
%     fnm = [resdir2 sessionID '_PSYCHPLOT.jpg'];   % save
%     saveas(H,fnm)
%     fnm = [resdir2 sessionID '_PSYCHPLOT.fig'];
%     saveas(H,fnm)
%     
%     H = auditory_gonogo_psychplot3(animalID,sessionID);
%     maximize_figure(H)
%     fnm = [resdir2 sessionID '_PSYCHPLOT2.jpg'];   % save
%     saveas(H,fnm)
%     fnm = [resdir2 sessionID '_PSYCHPLOT2.fig'];
%     saveas(H,fnm)
% end