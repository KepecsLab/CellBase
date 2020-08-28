function   [cellid, session_separator] = fname2cellid(fname,varargin)
%FNAME2CELLID    Convert filenames to cell IDs.
%   CELLID = FNAME2CELLID(FNAME) converts valid filenames into cellids
%   or returns empty if it fails.
%
%   CELLID = fname2cellid(FNAME,VARARGIN) uses VARARGIN{1} to look for
%   cellbase preferences instead of retrieving them using
%   getpref/getcbpref.
%
%   [CELLID, SESSION_SEPARATOR] = fname2cellid(...) additionally returns
%                          dtermined session separator
%
%   Valid filenames
%   (1) start with the default path or only include 'rat\session\unit.mat'
%   (2) unit 1 of tetrode 1 is called TT1_1.mat 
%   (3) session name can only contain '.' or '_' characters but not both
%   (4) should be consistent across the database
%
%   See alse FINDALLCELLS.

%   ZFM additions:
%   - Uses the preference 'cell_pattern' to store a search term that selects
%   files corresponding to units. The default is 'TT*.mat'. E.g.
%
%   setcbpref('cell_pattern','TT*.mat')
%
%   - Uses the preference 'group' to store an 'cut' subdirectory
%   under the session directory. E.g.
%
%   setcbpref('group','Mike');

% VARARGIN can either be empty or passing cellbase preference to avoid
% redunant getpref/getcbpref calls (TO).

%

%   Edit log: BH 3/21/11, TO 05/2018

% Get cellbase preferences
if isempty(varargin) %have not been passed --> get them directly
cellbase_path  = getpref('cellbase','datapath');
cell_pattern = getcbpref('Spikes_cell_pattern');
try
    getcbpref('group')
    group_param=true;
catch
    group_param=false;
end
else %prefs have been passed --> get them from varargin{1}
    cellbase_path = varargin{1}.datapath;
    cell_pattern = varargin{1}.Spikes_cell_pattern;
    if isfield(varargin{1},'group')
        group_param=true;
    else
        group_param=false;
    end
end
session_separator=''; %default

% Strip datpath from file
fn = char(strrep(fname,cellbase_path,''));
fs = filesep;

% Parse the filename (ratname\sessionname\analysisdir\)
[ratname,remain]  = strtok(fn,fs);
[session,remain]  = strtok(remain(2:end),fs);

% Added options for specifying an analysis directory below the session
% directory
if group_param
    % Strip analysis path
    [ad_junk,remain] = strtok(remain,fs);
end

% Extract tetrode unit
[tetrodeunit,ext] = strtok(remain(2:end),'.');
tu =  sscanf(tetrodeunit,[cell_pattern '%d_%d']);
pos_u = strfind(session,'_');
pos_p = strfind(session,'.');

% Control output
if ~strcmp(ext,'.mat')
    %disp('FNAM2CELLID: Not a matlab data file.');
    cellid = 0;
    return
end

if isempty(ratname) || isempty(session) || isempty(tu)
    strr = sprintf('FNAME2CELLID: Filename %s could not be parsed correctly.',fname);
   warning(strr) 
   cellid = 0;
   return
elseif ~isempty(pos_u) || ~isempty(pos_p)
    strr = sprintf('FNAME2CELLID: Filename %s could not be parsed correctly.',fname);
    warning(strr);
    cellid = 0;
    return    
elseif ~isempty(pos_u)   % there were underscores in the sessions
    session = strrep(session,'_','.');% replace them with .'s
    session_separator = '_'; 
%     setcbpref('session_separator','_');  % note this as a preference
elseif ~isempty(pos_p)
    session_separator = '.'; 
%     setcbpref('session_separator','.');
else
    % there is no separator, which is OK
    session_separator = ''; 
%     setcbpref('session_separator','');
end

cellid = sprintf('%s_%s_%d.%d',ratname,session,tu(1),tu(2));
%disp(cellid)