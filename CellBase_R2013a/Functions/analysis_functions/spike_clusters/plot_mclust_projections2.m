function HS = plot_mclust_projections2(cellid,varargin)
%PLOT_MCLUST_PROJECTIONS2   Plot feature data for clusters.
%   PLOT_MCLUST_PROJECTIONS2(CELLID) plots feature data for a tagged cell
%   (CELLID) and it's tetrode pairs in all projections. Tagged cluster is
%   shown in orange and light-evoked spikes are overlayed in blue. Time
%   period for light-evoked activity is selected automatically (see
%   FINDSTIMPERIOD). Only spikes from the beginning of the first to the end
%   of the last stimulation protocol are included. These default behaviors
%   can be modified using the following optional input arguments
%   (parameter, value pairs, with default values):
%       'stim_period', [] - start and end of stimulation period after each
%           light pulse
%       'feature_names', 'Energy' - features for which feature data are
%           plotted; character or cell array
%       'marker', '+' - marker for the scatter plots
%       'marker_size', 2 - marker size for the scatter plots
%       'usefastplot', true - use fast plotting method (downsample points
%           to plot only one point per pixel; appears the same); faster,
%           but zoom is not implemented (for saving in image formats, e.g.
%           pdf or jpg); if false, full data is plotted (for viewing or
%           saving fig format)
%       'stimonly', true - only spikes from the beginning of the first to
%           the end of the last stimulation protocol are selected for
%           plotting; if false, all spikes are included
%       'plotlightspikes', true - if true, light-evoked spikes are 
%           superimposed
%       'fighandle', NaN - optional figure handle for feature plots
%       'axhandle', NaN - optional axes handles for individual feature
%           plots
%       'plot', true - determines if features are plotted
%       'plotbest', false - only plots best feature combination based on a
%          simple metrix
%
%   HS = PLOT_MCLUST_PROJECTIONS2(CELLID,...) returns the handles of the
%   figures in HS struct. The fields of HS are named according to the
%   features and channels: HS.(XFeatureXChannel_YFeatureYChannel), e.g.
%   HS.Amplitude1_Energy4.
%
%   Example:
%   plot_mclust_projections2(cellid,'feature_names',{'Amplitude','Energy'},...
%        'stim_period',[0.002 0.004]);
%
%   See also PLOTWAVEFORMS.

%   Sachin Ranade & Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   08-May-2012

%   Edit log: SPR 12/28/11; BH 5/8/12; TO 12/2017

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'stim_period',[],@isnumeric)   % start and end point for the interval of light-evoked spikes
addParamValue(prs,'feature_names',{'Energy'},@(s)iscell(s)|ischar(s))  % names of MClust features to use
addParamValue(prs,'marker','+')  % marker for the plots
addParamValue(prs,'marker_size',2,@isnumeric)  % marker size for the plots
addParamValue(prs,'usefastplot',true,@(s)islogical(s)|ismember(s,[0 1]))  % fast plotting (no zoom)
addParamValue(prs,'stimonly',true,@(s)islogical(s)|ismember(s,[0 1]))  % restrict to period between first and last light pulse
addParamValue(prs,'plotlightspikes',true,@(s)islogical(s)|ismember(s,[0 1]))  % plot light-evoked spikes
addParamValue(prs,'fighandle',NaN,@(s) isnumeric(s)|isgraphics(s,'figure'))  % figure handle for feature plots
addParamValue(prs,'axhandle',NaN,@(s)isgraphics(s,'axes')|isnumeric(s))  % axes handles for individual feature plots
addParamValue(prs,'plot',true,@(s)islogical(s)|ismember(s,[0 1]))  % axes handles for individual feature plots
addParamValue(prs,'plotbest',false,@(s)islogical(s)|ismember(s,[0 1]))  % axes handles for individual feature plots
parse(prs,cellid,varargin{:})
g = prs.Results;
if ischar(g.feature_names)
    g.feature_names = {g.feature_names};
end
if g.plotbest
    g.plot=false;
end

% Load stim events to get Pulse onset times.
ST = loadcb(cellid,'Stimevents');
pon = ST.PulseOn(~isnan(ST.PulseOn));

% Load spikes from Ntt file.
Nttfn = cellid2fnames(cellid,'Ntt');
loadingEngine = getcbpref('TrodeLoadingEngine');
[all_spikes,~] = feval(loadingEngine,Nttfn);
blocks = [find(diff(pon)>1600),length(pon)];
if length(blocks)==1
    val_spk=all_spikes >= pon(1) & all_spikes <= pon(blocks(end));
else
    val_spk=all_spikes >= pon(1) & all_spikes <= pon(blocks(1));
for b = 1:length(blocks)-1
    val_spk = val_spk | (all_spikes>=pon(blocks(b)+1) & all_spikes<=pon(blocks(b+1)));
end
end
val_spk_i = find(val_spk); % consider spikes only within the stimulation protocol/blocks to account for drift
nspk = length(all_spikes);
TIMEFACTOR = getcbpref('Spikes_timefactor');    % scaling factor to convert spike times into seconds
spk = loadcb(cellid,'Spikes'); spk = spk*TIMEFACTOR;
[junk,junk2,tagged_cell_inx] = intersect(spk,all_spikes);  %#ok<*ASGLU> % get indices for the cell
if ~isequal(junk,spk)  % check if all files have appropriate time stamps
    error('plot_mclust_projections:SpikeTimeMismatch','Mismatch between saved spike times and Ntt time stamps.')
end
if g.stimonly   % restrict to stimulation epoch
    tagged_cell_inx = intersect(tagged_cell_inx,val_spk_i);
end
tagged_cell_i = zeros(nspk,1);
tagged_cell_i(tagged_cell_inx) = 1;

% Cells on same tetrode including cellid.
[NumCell,tetpartners] = tetrodepairs(cellid);
tag_cell = strmatch(cellid,tetpartners);

% Spikes from each cell have an index. Spikes from noise have index 0.
cell_i = zeros(nspk,1);
for iCell = 1:NumCell
    spk = loadcb(tetpartners(iCell),'Spikes');spk=spk*10^-4; % load spike times.
    [junk,junk2,cell_inx] = intersect(spk,all_spikes); % get indices for the cell
    if g.stimonly   % restrict to stimulation epoch
        cell_inx = intersect(cell_inx,val_spk_i);
    end
    cell_i(cell_inx) = iCell;
end
%restrict noise spikes to stimulation epoch
if g.stimonly   % restrict to stimulation epoch
    all_i = 1:nspk;
    invalid_spk_i = setxor(all_i,val_spk_i);
    cell_i(invalid_spk_i) = -1;
end

% Load feature data for tetrode.
[r,s,t] = cellid2tags(cellid);
for k = 1:length(g.feature_names)
    prop = [g.feature_names{k} '.fd'];
    propfn = [getcbpref('Spikes_cell_pattern') num2str(t) '_' prop];
    sessionpath = cellid2fnames(cellid,'sess');
    propfn_path = [sessionpath filesep 'FD'];   % where the feature file can be found
    if ~isdir(propfn_path)
        propfn_path = sessionpath;
    end
    propfn_path = fullfile(propfn_path,propfn);
    wf_prop = load(propfn_path,'-mat');
    FeatureData(k,:,:) = wf_prop.FeatureData; %#ok<AGROW>
end

% Light-evoked spikes
if g.plotlightspikes
    
    % Latency of stimulated spikes
    if isempty(g.stim_period)
        [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
        if isnan(lim1) || isnan(lim2)
            lim1 = 0.001;   % no activation detected
            lim2 = 0.004;
        end
    else
        lim1 = g.stim_period(1);
        lim2 = g.stim_period(2);
    end
    
    % Evoked spikes
    tsegs_evoked = rel2abstimes(cellid,[lim1 lim2],'stim','PulseOn');   % convert period to epochs relative to pulses
    selts_evoked = extractSegSpikes(all_spikes,tsegs_evoked);   % find putative stimualated spikes
    [junk,junk2,evoked_cell_inx] = intersect(selts_evoked,all_spikes); % get indices for light-evoked spikes
end

% Plot
%plot features (only if g.plot=true)
if g.plot
    HS = plot_features(g,NumCell,cell_i,tagged_cell_i,evoked_cell_inx,FeatureData,'all');
else
    HS = plot_features(g,NumCell,cell_i,tagged_cell_i,evoked_cell_inx,FeatureData,'none');%get metric
end

if g.plotbest
    if  ~isgraphics(g.axhandle(1),'axes')
        ax=axes;
    else
        ax = g.axhandle(1);
        axes(ax)
    end
    HS = plot_features(g,NumCell,cell_i,tagged_cell_i,evoked_cell_inx,FeatureData,HS.best_feature_idx);
    
end
end
    

% -------------------------------------------------------------------------
function HS = plot_features(g,NumCell,cell_i,tagged_cell_i,evoked_cell_inx,FeatureData,whichcmb)
cmp = hsv(NumCell) / 4 + 0.75;
pcmb = allcomb(1:size(FeatureData,1),1:size(FeatureData,3));
cmb = flipud(combnk(pcmb(:,1)*10+pcmb(:,2),2));
NumComb = size(cmb,1);
plotc=g.plot|g.plotbest;
if strcmpi(whichcmb,'all')
    Comb = 1:NumComb;
elseif strcmp(whichcmb,'none') || isnan(whichcmb)
    Comb = 1:NumComb;
    plotc=false;
else
   Comb = whichcmb;
   NumComb = length(Comb);
   g.axhandle(Comb) = g.axhandle(1:NumComb);
end

NCols = ceil(NumComb/3);
NRows = ceil(NumComb/NCols);
ax=zeros(1,NumComb);
dist_metric=zeros(1,NumComb);
%open figure
if plotc && ~isnumeric(g.fighandle) && isgraphics(g.fighandle,'figure')
    %valid fig handle --> plot all projections in figure
    HS.feature_plot=figure(g.fighandle);
    plot_to_figs=false;
elseif plotc && any(~isnumeric(g.axhandle)) && any(isgraphics(g.axhandle,'axes') )
    plot_to_figs=false;
else%no valid fig or ax handle --> plot projections to individual figures
    plot_to_figs=true;
end
k_i=0;
for k = Comb
    k_i=k_i+1;
    fst = [floor(cmb(k,1)/10) mod(cmb(k,1),10)];
    scnd = [floor(cmb(k,2)/10) mod(cmb(k,2),10)];
    xdata = squeeze(FeatureData(fst(1),cell_i==0,fst(2)));
    ydata = squeeze(FeatureData(scnd(1),cell_i==0,scnd(2)));
    
    % select plot
    titlestr = [g.feature_names{fst(1)} num2str(fst(2)) '_'...
        g.feature_names{scnd(1)} num2str(scnd(2))];
    if plotc
        if ~plot_to_figs && (length(g.axhandle)<k || ~isgraphics(g.axhandle(k),'axes'))
            ax(k)=subplot(NRows,NCols,k_i);
        elseif ~plot_to_figs
            ax(k) = g.axhandle(k);
            axes(ax(k))
        else %plot_to_figs
            HS.feature_plot(k)=figure;
            ax(k) = axes;
        end
    end
    hold on
    axis([min(xdata) max(xdata) min(ydata) max(ydata)])
    xlabel([g.feature_names{fst(1)} ': ' num2str(fst(2))])
    ylabel([g.feature_names{scnd(1)} ': ' num2str(scnd(2))])
    
    % Plot noise spikes
    if g.usefastplot && plotc
        fastplot(xdata,ydata,[0.8 0.8 0.8],g.marker,g.marker_size);
    elseif plotc
        slowplot(xdata,ydata,[0.8 0.8 0.8],g.marker,g.marker_size);
    end
    
    % Plot all clusters
    for iC = 1:NumCell
        xdatai = squeeze(FeatureData(fst(1),cell_i==iC,fst(2)));
        ydatai = squeeze(FeatureData(scnd(1),cell_i==iC,scnd(2)));
        if g.usefastplot && plotc
            fastplot(xdatai,ydatai,cmp(iC,:),g.marker,g.marker_size);
        elseif plotc
            slowplot(xdatai,ydatai,cmp(iC,:),g.marker,g.marker_size);
        end
    end
    
    % Plot tagged cluster
    xdatai = squeeze(FeatureData(fst(1),tagged_cell_i==1,fst(2)));
    ydatai = squeeze(FeatureData(scnd(1),tagged_cell_i==1,scnd(2)));
    if g.usefastplot && plotc
        fastplot(xdatai,ydatai,[255 204 0]/255,g.marker,g.marker_size);
    elseif plotc
        slowplot(xdatai,ydatai,[255 204 0]/255,g.marker,g.marker_size);
    end
    
    %separability metric
    c_noise = [xdata',ydata'];
    c_noise = (cov(c_noise)^-1)*c_noise';
    c_noise_m = median(c_noise,2);
    c_tagged = [xdatai',ydatai'];
    c_tagged = (cov(c_tagged)^-1)*c_tagged';
    c_tagged_m = median(c_tagged,2);
    dist_metric(k) =sum(abs(c_noise_m-c_tagged_m));
    
    % Plot light-evoked spikes
    if g.plotlightspikes
        xdatai = squeeze(FeatureData(fst(1),evoked_cell_inx,fst(2)));
        ydatai = squeeze(FeatureData(scnd(1),evoked_cell_inx,scnd(2)));
        if g.usefastplot && plotc
            fastplot(xdatai,ydatai,[0 153 255]/255,'.',g.marker_size+5);
        elseif plotc
            slowplot(xdatai,ydatai,[0 153 255]/255,'.',g.marker_size+5);
        end
    end
end
HS.dist_metric = dist_metric;
[~,best_i] = max(dist_metric);
fst = [floor(cmb(best_i,1)/10) mod(cmb(best_i,1),10)];
    scnd = [floor(cmb(best_i,2)/10) mod(cmb(best_i,2),10)];
best_name = [g.feature_names{fst(1)} num2str(fst(2)) '_'...
        g.feature_names{scnd(1)} num2str(scnd(2))];
HS.best_feature = best_name;
HS.best_feature_idx = best_i;
if plotc && g.plot && length(ax)>=k && isgraphics(ax(k),'axes')
    set(ax(best_i),'Color',[1,.95,.95]);
end
end

function C = allcomb(A,B)

% Convert to columns
A = A(:);
B = B(:);

% Combinations
as = arrayfun(@(k)horzcat(repmat(A(k),length(B),1),B),1:length(A),'UniformOutput',false);
C = cell2mat(as');
end

% -------------------------------------------------------------------------
function h = fastplot(x,y,clr,mrk,mrks)

% Reduce number of points
old_units = get(gca,'Units');
set(gca,'Units','pixels')
pos = get(gca,'Position');
xpixels = pos(3) + 1;
ypixels = pos(4) + 1;

xl = xlim;
mnx = xl(1);
mxx = xl(2);
yl = ylim;
mny = yl(1);
mxy = yl(2);
x2 = round((x-mnx)/(mxx-mnx)*(xpixels-1)) + 1;
y2 = round((y-mny)/(mxy-mny)*(ypixels-1)) + 1;
u = unique(x2*100000+y2);
y3 = mod(u,100000);
x3 = (u - y3) / 100000;
x4 = (x3 / xpixels) * (mxx - mnx) + mnx;
y4 = (y3 / ypixels) * (mxy - mny) + mny;

% Plot
h = plot(x4,y4,'.','Color',clr,'Marker',mrk,'MarkerSize',mrks);

% Restore axis units property
set(gca,'Unit',old_units)
end

% -------------------------------------------------------------------------
function h = slowplot(x,y,clr,mrk,mrks)

h = plot(x,y,'.','Color',clr,'Marker',mrk,'MarkerSize',mrks);
end