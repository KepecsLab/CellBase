function error_list = taggedsummaryTO(I,issave)
%TAGGEDSUMMARY   Summary properties of putative tagged neurons.
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

%   Edit log: BH 5/10/12, 4/19/13; TO 12/2017

% Pass the control to the user in case of error
 
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
        try
            HPO = Hindex(findcellpos(cellid)); % already calculated by TAGGING, so pass it on
        catch %if not, call tagging
            tagging(cellid);
            Hindex = getvalue('Hindex');
            HPO = Hindex(findcellpos(cellid));
        end    
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
% ST = loadcb(cellid,'STIMSPIKES');
% if isequal(findcellstr(ST.events(:,1),'PulseOnS'),0)
%     [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
%     prealignSpikes(cellid,'events',...
%         {'PulseOnS' 'PulseOn' 'PulseOn' [lim1 lim2]},...
%         'epochs',[],'filetype','stim','ifsave',1,'ifappend',1)
% else
%     [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
%     prealignSpikes(cellid,'events',...
%         {'PulseOnS' 'PulseOn' 'PulseOn' [lim1 lim2]},...
%         'epochs',[],'filetype','stim','ifsave',1,'writing_behavior','replace')
% end

SE = loadcb(cellid,'StimEvents');


% Efficiency, latency and jitter for 'PulseOn'
[E_pulseon L_pulseon J_pulseon B_pulseon M_pulseon A1 A2] = ...
    reliability_latency_jitter(cellid,'event','PulseOn','fighandle',NaN,'axhandle',NaN);

% H-index for 'PulseOn'
Hindex_pulseon = HPO;

%Plot
fighandle = figure('Position',[ 488.2000  165.8000  883.2000  596.0000],'NumberTitle','off','MenuBar','none','Color',[1,1,1],'Name','Tagging summary');

% BurstOn and PulseOn PSTH
psth_ax = subplot(4,4,[1,2,5,6]);
psth_ax_back = plot_raster_psth(cellid,'BurstOn',false,'PulseOn',true,'PulseFig',fighandle.Number,'PulseAx',psth_ax);
psth_ax = psth_ax_back.HA_PulseOn{end}; %last one is for psth


% Plot light-evoked and spont. spikes
wf_ax(1)=subplot(4,4,9);
wf_ax(2)=subplot(4,4,10);
wf_ax(3)=subplot(4,4,13);
wf_ax(4)=subplot(4,4,14);
HWF = plotwaveforms(cellid,'correlation',true,'maxnum',5000,'fighandle',[NaN,NaN,fighandle.Number],'axhandle',{[],[],wf_ax},'evoked',false,'spont',false);
for i =1:length(wf_ax)
     wf_ax(i).XAxis.Visible='off';
     wf_ax(i).YAxis.Visible='off';
end



% Plot MClust projections
proj = subplot(4,4,[11,12,15,16]);
plot_mclust_projections2(cellid,'feature_names',{'Peak','Energy'},...
    'stim_period',[A1 A2],'stimonly',true,'usefastplot',true,'plot',true,'fighandle',fighandle.Number,'axhandle',proj,'plotbest',true);
MinMaxLim(proj,0.995)
proj.XAxis.TickValues=[];
proj.YAxis.TickValues=[];

%write statistics
%plot
yl = psth_ax.YLim;
if ~isnan(A1) && ~isnan(A2)
    period = fill(psth_ax,[A1,A2,A2,A1],[yl(1),yl(1),yl(2),yl(2)],[1,.9,.9]);
    period.EdgeColor='none';
    psth_ax.Legend.String = psth_ax.Legend.String(1:end-1);
    psth_ax.Children = psth_ax.Children([2:end,1]);
end
l_line = line(psth_ax,[L_pulseon,L_pulseon],yl','Color',[.8,0,0],'LineWidth',1);
l_line_j = line(psth_ax,[L_pulseon-J_pulseon,L_pulseon+J_pulseon],[yl(2),yl(2)],'Color',[.8,0,0],'LineWidth',.5);
psth_ax.Legend.String = psth_ax.Legend.String(1:end-2);

%write
[lID_amp, lLr_amp, Pref ,Pref2] = lightcluster(cellid,'feature_names',{'Peak','Energy'},...
    'stim_period',[A1 A2]);
[lID_PC, lLr_PC] = lightcluster(cellid,'feature_names',{'Peak','WavePC1'},...
    'stim_period',[A1 A2]);

leftalign = 0.6;
rightalign = 0.7;
leftwidth=0.2;
rightwidth = 0.2;
height=0.3;

str_head = {'Cluster quality measures';};
str_left = {'Latency: ';...
    'Jitter: ';...
    'Reliability: ';...
    'H index: ';...
    ' ';...
    'ID (A, E): ';...
    'L-ratio (A, E): ';...
    'ID (PC, E): ';...
    'L-ratio (PC, E): ';...
    'Preference 1: ' ;...
    'Preference 2: ';...
    'R-Waveform: ';...
    };
str_right =  {   [num2str(L_pulseon*1000) ' ms'];...
    [num2str(J_pulseon*1000) ' ms'];...
    [ num2str(E_pulseon*100) ' %'];...
    [ num2str(Hindex_pulseon) ];...    
    ' ';...
    [ num2str(lID_amp)];...
    [ num2str(lLr_amp)];...
    [ num2str(lID_PC)];...
    [ num2str(lLr_PC)];...
    [ num2str(Pref)];...
    [ num2str(Pref2)];
    [ num2str(HWF.R)]};

uicontrol(fighandle,'Style','text','Unit','normalized','Position',...
    [leftalign, 0.85, 0.2, 0.05],'FontSize',12,'HorizontalAlignment',...
    'left','String',str_head,'BackgroundColor',[1 1 1],'FontWeight','bold');

uicontrol(fighandle,'Style','text','Unit','normalized','Position',...
    [leftalign,0.55 ,leftwidth, height],'FontSize',12,'HorizontalAlignment',...
    'left','String',str_left,'BackgroundColor',[1 1 1]);

uicontrol(fighandle,'Style','text','Unit','normalized','Position',...
    [rightalign, 0.55, rightwidth, height],'FontSize',12,'HorizontalAlignment',...
    'left','String',str_right,'BackgroundColor',[1 1 1]);

set(findall(fighandle,'-property','FontName'),'FontName','Helvetica')

% Save
if issave
    if ~isdir(resdir),mkdir(resdir),end
%     save(fullfile(resdir,strcat('TAGGEDPROP_', regexprep(cellid,'\.','_'), '.mat')),...
%         'Hindex_burston','D_KL_burston',...
%         'E_pulseon','L_pulseon','J_pulseon','B_pulseon','M_pulseon','A1','A2',...
%         'E_burston','L_burston','J_burston','B_burston','M_burston',...
%         'lLr_amp','lID_amp','lLr_PC','lID_PC','Pref','Pref2')
    
    % Write figures to pdf
    pdfname = fullfile(resdir,'TAGGEDSUMMARY.pdf');
    writefigs(fighandle,pdfname)
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
        if ishandle(h)
            export_fig(h,'-dpdf',pdfname,'-painters','-append');  % write to pdf
        end
    end
else
    export_fig(H,'-dpdf',pdfname,'-painters','-append');  % write to pdf
end