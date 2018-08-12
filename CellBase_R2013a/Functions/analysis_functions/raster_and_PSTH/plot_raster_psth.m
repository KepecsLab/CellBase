function HS = plot_raster_psth(cellid,varargin)
%PLOT_RASTER_PSTH   Plot raster plots and PSTHs.
%   HS = PLOT_RASTER_PSTH(CELLID) plots raster plots and PSTHs aligend to
%   'BurstOn' and/or 'PulseOn' events. Figure handles are returned in HS.
%   Which rasters/PSTHs are plotted is controlled by the following optional
%   input argument parameter-value pairs (with default values):
%       'BurstOn', true - if true, plot 'BurstOn' raster and PSTH
%       'PulseOn', false - if true, plot 'PulseOn' raster and PSTH
%       'BurstFig', NaN - figure handle to plot BurstOn
%       'BurstAx', NaN - axes handle to plot BurstOn
%       'PulseFig', NaN - figure handle to plot PulseOn
%       'PulseAx', NaN - axes handle to plot PulkseOn
%
%   Output (fields of HS struct):
%       H_BurstOn - figure handle for 'BurstOn' raster and PSTH
%       H_PulseOn - figure handle for 'PulseOn' raster and PSTH
%
%   See also VIEWCELL2B and PLOT_RASTER2A.

%   Edit log: BH 5/9/12; TO 12/2017

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'BurstOn',true,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying 'BurstOn' rasters and PSTHs
addParamValue(prs,'PulseOn',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying 'PulseOn' rasters and PSTHs
addParamValue(prs,'BurstFig',NaN,@isnumeric)   % control plotting
addParamValue(prs,'PulseFig',NaN,@isnumeric)   % control plotting
addParamValue(prs,'BurstAx',NaN,@(s)isnumeric(s)|isgraphics(s,'axes'))   % control plotting
addParamValue(prs,'PulseAx',NaN,@(s)isnumeric(s)|isgraphics(s,'axes'))   % control plotting
parse(prs,cellid,varargin{:})
g = prs.Results;

% Set input parameters for 'viewcell2b'
SEvent = 'BurstOn';
winBurst = [-0.2 2];
winPulse = [-0.025,0.05];
partsPO = 'all';
%partsBO = '#BurstNPulse';
partsBO='all';
dt = 0.001;
sigma = 0.001;
PSTHstd = 'on';
% ShEvent = {'PulseOn','BurstOff'};
ShEvent = {'PulseOn'}; % FS MOD, Kludge 
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);

% Plot raster plot and PSTH for 'BurstOn'
if g.BurstOn
    if ~isempty(g.BurstFig) && ~isnan(g.BurstFig(1))
        HS.H_BurstOn = figure(g.BurstFig(1));
    else
        HS.H_BurstOn = figure();
    end
    if ~isempty(g.BurstAx) && isgraphics(g.BurstAx(1),'axes')
        HS.HA_BurstOn = g.BurstAx(1);
    else
        HS.HA_BurstOn = axes;
    end     
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    HS.HA_BurstOn = viewcell2b(cellid,'TriggerName','BurstOn','SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',ShEvColors,...
        'FigureNum',HS.HA_BurstOn,'eventtype','stim','window',winBurst,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',partsBO,...
        'EventMarkerWidth',0,'PlotZeroLine','off');
    HS.HA_BurstOn = [];
%     for i =1:length(axh) %hack to return axes objects (since M2014f everything is in objects, best work with that)
%         axes(axh(i));
%         HS.HA_BurstOn{i} = gca;
%     end    
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end

% Plot raster plot and PSTH for 'PulseOn'
if g.PulseOn
    if ~isempty(g.PulseFig) && ~isnan(g.PulseFig(1))
        HS.H_PulseOn = figure(g.PulseFig);
    else
        HS.H_PulseOn = figure('Position',[97 163 1533 815]);
    end
    if ~isempty(g.PulseAx) && isgraphics(g.PulseAx(1),'axes')
        HS.HA_PulseOn = g.PulseAx(1);
    else
        HS.HA_PulseOn = axes;
    end    
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    axh = viewcell2b(cellid,'TriggerName','PulseOn','SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',ShEvColors,...
        'FigureNum',HS.HA_PulseOn,'eventtype','stim','window',winPulse,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',partsPO,...
        'EventMarkerWidth',0,'PlotZeroLine','off','Num2Plot',1000);
    HS.HA_PulseOn = [];
%     for i =1:length(axh) %hack to return axes objects (since M2014f everything is in objects, best work with that)
%         axes(axh(i));
%         HS.HA_PulseOn{i} = gca;
%     end
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end