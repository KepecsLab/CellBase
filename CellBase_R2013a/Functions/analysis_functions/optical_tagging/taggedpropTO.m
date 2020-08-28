function error_list = taggedpropTO(I,issave)
%TAGGEDPROP   Properties of putative tagged neurons.
%   TAGGEDPROP(I,ISSAVE) performs analysis aiding the definite decision
%   about tagging. Input parameters: 
%       I - list of cell IDs or index set to CELLIDLIST (see CellBase
%           documentation); if empty or not specified, putative tagged
%           cells are selected (ID>20, L-ratio<0.15, H index<0.01, R>0.9;
%           see TAGGING, LRATIO, NBISSTIM and SPIKESHAPECORR for
%           details on these measures)
%       ISSAVE - controls saving
%
%   The following analyses are performed (output variables saved in mat
%   files and figures in pdf):
%       Reliability, latency and jitter of evoked spikes after 'BurstOn',
%           'PulseOn' and frequency-restricted 'PulseOn' events; see
%           RELIABILITY_LATENCY_JITTER for details
%       H index for 'BurstOn', 'PulseOn' and frequency-restricted 'PulseOn'
%           events; see TAGGING_INDEX for details
%       Waveforms of spontaneous and light-evoked spikes; see PLOTWAVEFORMS
%           for details
%       All projections of feature data in the Energy-Amplitue space with
%           the putative tagged cells shown in orange and the light-evoked 
%           spikes overlayed in blue; see PLOT_MCLUST_PROJECTIONS2 for
%           details
%       Cluster quality indices restricted to light-evoked spikes in the
%           Energy-Amplitude as well as in the Energy-WavePC1 space; see
%           LIGHTCLUSTER and LRATIO for details
%       Raster plot and PSTH aligned to 'BurstOn' and 'PulseOn' events
%           (only 1000 pulses shown in the latter); see PLOT_RASTER_PSTH,
%           VIEWCELL2B and PLOT_RASTER2A for details.
%
%   ERROR_LIST = TAGGEDPROP(I,ISSAVE) returns a structure with all caught
%   errors. ERROR_LIST includes cell IDs, error messages and the captured
%   exception variables. A time stamped ERROR_LIST is also assigned in base
%   workspace and saved to the results directory automatically.
%
%   See also TAGGING, LRATIO, SPIKESHAPECORR, RELIABILITY_LATENCY_JITTER,
%   TAGGING_INDEX, PLOTWAVEFORMS, PLOT_MCLUST_PROJECTIONS2, LIGHTCLUSTER,
%   PLOT_RASTER_PSTH, VIEWCELL2B and PLOT_RASTER2A.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   05-Oct-2012

%   Edit log: BH 5/10/12, 4/19/13

% Pass the control to the user in case of error
dbstop if error
 
% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = true;
end
if nargin < 1
    I = [];
end
 
% Directories
cbpath = getpref('cellbase','datapath');
resdir = fullfile(cbpath,'taggedprop');
if ~isdir(resdir)
    mkdir(resdir)
end
 
% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

% Find putative tagged cells
if isempty(I)
    Lratio = getvalue('Lr_PC');
    ID = getvalue('ID_PC');
    Hindex = getvalue('Hindex');
    R = getvalue('R');
    ptinx = ID > 20 & Lratio < 0.15 & Hindex < 0.01 & R > 0.9;
    I = find(ptinx);
    putative_tagged = CELLIDLIST(I);
else
    Hindex = getvalue('Hindex');
    if isnumeric(I)
        putative_tagged = CELLIDLIST(I);   % index set to CELLIDLIST
    elseif ischar(I)
        putative_tagged = {I};   % only one cellID passed
    elseif iscellstr(I)
        putative_tagged = I;   % list of cell IDs
    else
        error('taggedprop:inputArg','Unsupported format for cell IDs.')
    end
end

% All the analysis promised in the help
NumCells = length(putative_tagged);
error_list = struct('cellid',{},'message',{},'exception',{});  % keep a list of caught errors 
errinx = 0;
for k = 1:NumCells
    cellid = putative_tagged{k};
    try
        HPO = Hindex(findcellpos(cellid)); % already calculated by TAGGING, so pass it on
        main(cellid,issave,resdir,HPO)  % every analysis
    catch ME
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)
        errinx = errinx + 1;  % error counter
        error_list(errinx).cellid = cellid;   % record problematic cell ID
        error_list(errinx).message = ME.message;   % error message
        error_list(errinx).exception = ME;   % exception structure
    end
end

% Create time-stamped error list in base workspace
dsr = datestr(now);  % date stamp
dsr = regexprep(dsr,':','_');
list_name = ['error_list_' dsr];  % variable name
list_name = regexprep(list_name,'-','_');
list_name = regexprep(list_name,' ','_');
assignin('base',list_name,error_list)   % assign error list in base workspace
error_fname = fullfile(resdir,[list_name '.mat']);   % file name
save(error_fname,'error_list')   % save error list

% -------------------------------------------------------------------------
function main(cellid,issave,resdir,HPO)

% Add 'PulseOn' event if missing
ST = loadcb(cellid,'STIMSPIKES');
if isequal(findcellstr(ST.events(:,1),'PulseOnS'),0)
    [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
    prealignSpikes(cellid,'events',...
        {'PulseOnS' 'PulseOn' 'PulseOn' [lim1 lim2]},...
        'epochs',[],'filetype','stim','ifsave',1,'ifappend',1)
else
    [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
    prealignSpikes(cellid,'events',...
        {'PulseOnS' 'PulseOn' 'PulseOn' [lim1 lim2]},...
        'epochs',[],'filetype','stim','ifsave',1,'writing_behavior','replace')
end

SE = loadcb(cellid,'StimEvents');


% Efficiency, latency and jitter for 'PulseOn'
[E_pulseon L_pulseon J_pulseon B_pulseon M_pulseon A1 A2] = ...
    reliability_latency_jitter(cellid,'event','PulseOn','fighandle',NaN,'axhandle',NaN);

% H-index for 'PulseOn'
Hindex_pulseon = HPO;

% Efficiency, latency and jitter for 'BurstOn'
[E_burston L_burston J_burston B_burston M_burston] = ...
    reliability_latency_jitter(cellid,'event','BurstOn');

% H-index for 'BurstOn'
[Hindex_burston D_KL_burston] = tagging_index(cellid,'event','BurstOn');  %#ok<NASGU> % H-index, D_KL



% Plot efficiency
HE = figure('Position',[624 110 900 868]);
S = set_subplots(4,1,0.065,0.065);
axes(S(1)) %#ok<*MAXES>
L1 = line(xlim,[E_pulseon E_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
L2 = line(xlim,[E_burston E_burston],'Color',[255 204 0]/255,'LineWidth',2);
legend([L1 L2],{'PulseOn' 'BurstOn'},'Location','EastOutside')
title('Efficiency')

% Plot latency and jitter
axes(S(2))
hold on
L1 = line(xlim,[L_pulseon L_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
line(xlim,[L_pulseon+J_pulseon L_pulseon+J_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle','--');
line(xlim,[L_pulseon-J_pulseon L_pulseon-J_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle','--');
L2 = line(xlim,[L_burston L_burston],'Color',[255 204 0]/255,'LineWidth',2);
line(xlim,[L_burston+J_burston L_burston+J_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle','--');
line(xlim,[L_burston-J_burston L_burston-J_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle','--');
legend([L1 L2],{'PulseOn' 'BurstOn'},'Location','EastOutside')
title('Latency')

% Plot evoked firing rate
axes(S(3))

L1 = line(xlim,[M_pulseon M_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
L2 = line(xlim,[M_burston M_burston],'Color',[255 204 0]/255,'LineWidth',2);
L3 = line(xlim,[B_pulseon B_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle','--');
L4 = line(xlim,[B_burston B_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle','--');
legend([L1 L2 L3 L4],{'PulseOn' 'BurstOn' 'Baseline' 'Baseline'},...
    'Location','EastOutside')
title('Evoked firing rate')

% Plot H-index
axes(S(4))

L1 = line(xlim,[Hindex_pulseon Hindex_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
L2 = line(xlim,[Hindex_burston Hindex_burston],'Color',[255 204 0]/255,'LineWidth',2);
line(xlim,[0.01 0.01],'Color','r','LineWidth',2)
legend([L1 L2],{'PulseOn' 'BurstOn'},'Location','EastOutside')
title('H index')

% Plot light-evoked and spont. spikes
HW = plotwaveforms(cellid,'correlation',true,'maxnum',5000,'fighandle',NaN,'axhandle',{});

% Plot MClust projections]
if HPO>0.1
    lightlim = [0.001,0.004];
else
    lightlim=[];%determined by peak
end
HM = plot_mclust_projections2(cellid,'feature_names',{'Peak','Energy'},...
    'stim_period',[A1 A2],'stimonly',true,'usefastplot',true,'plot',true,'fighandle',NaN,'axhandle',NaN,'plotbest',false,'stim_period',lightlim);

% Distance from light-evoked noise
[lID_amp lLr_amp Pref Pref2] = lightcluster(cellid,'feature_names',{'Peak','Energy'},...
    'stim_period',[A1 A2]);
[lID_PC lLr_PC] = lightcluster(cellid,'feature_names',{'Peak','WavePC1'},...
    'stim_period',[A1 A2]);
HD = figure;
axes;
str = {'Cluster quality measures for light-evoked spikes:';...
    ' ';...
    ['ID (Amplitude, Energy): ' num2str(lID_amp)];...
    ['L-ratio (Amplitude, Energy): ' num2str(lLr_amp)];...
    ['ID (WavePC1, Energy): ' num2str(lID_PC)];...
    ['L-ratio (WavePC1, Energy): ' num2str(lLr_PC)];...
    ['preference: ' num2str(Pref)];...
    ['preference: ' num2str(Pref2)]};
uicontrol('Style','text','Unit','normalized','Position',...
    [0.18 0.3 0.65 0.5],'FontSize',12,'HorizontalAlignment',...
    'left','String',str,'BackgroundColor',[1 1 1]);
axis off

% BurstOn and PulseOn PSTH
HR = plot_raster_psth(cellid,'BurstOn',false,'PulseOn',true,'PulseFig',NaN,'PulseAx',NaN);

% Save
if issave
    if ~isdir(resdir),mkdir(resdir),end
    save(fullfile(resdir,strcat('TAGGEDPROP_', regexprep(cellid,'\.','_'), '.mat')),...
        'Hindex_burston','D_KL_burston',...
        'E_pulseon','L_pulseon','J_pulseon','B_pulseon','M_pulseon','A1','A2',...
        'E_burston','L_burston','J_burston','B_burston','M_burston',...
        'lLr_amp','lID_amp','lLr_PC','lID_PC','Pref','Pref2')
    
    % Write figures to pdf
    pdfname = fullfile(resdir,['TAGGEDPROP_' regexprep(cellid,'\.','_') '.pdf']);
    writefigs(HR,pdfname)
    writefigs(HE,pdfname)
    writefigs(HW,pdfname)
    writefigs(HM,pdfname)
    writefigs(HD,pdfname)
    close all
end

% -------------------------------------------------------------------------
function writefigs(H,pdfname)
% Write figures to pdf
 
% Append to pdf
if isstruct(H)  % H is either a struct with figure handles, or a single fig. handle
    fls = fieldnames(H);
    for fs = 1:length(fls)
        h = H.(fls{fs});
        for i =1:length(h)
            if ishandle(h(i))
                export_fig(h(i),'-append',pdfname);  % write to pdf
            end
        end
    end
else
    export_fig(H,'-append',pdfname);  % write to pdf
end