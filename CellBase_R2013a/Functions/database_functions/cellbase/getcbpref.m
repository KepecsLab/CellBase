function val = getcbpref(varargin)
% GETCBPREF returns current cellbase's preference value(s).
% Cellbase-specific preferences are stored in cellbase file
% (getpref('cellbase','fname')). Note that cellbase-general preferences
% (cellbase name, data folder, cellbase file) are manipulated with
% getpref('cellbase')/setpref('cellbase',...).
%
% Usage
%
% val = getcbpref() - returns all cellbase preferences in struct VAL (including global preferenes).
% val = getcbpref(name) - returns cellbase preference with name NAME in VAL.
%
% Example
%
% val = getcbpref('TrialEvents_fname')
%
% See also: setcbpref(name,value), default_preferences()
%
% TO 05/2018

if nargin<2
    if nargin==1 && any(strcmp({'datapath','name','fname','filesep','cellbases'},varargin{1})) %global preferencces
        warning('Please use getpref to get global cellbase settings.');
        val = getpref('cellbase',varargin{1});     
    else %correct use for cellbase-specific preferences
        if exist(fullfile(getpref('cellbase','fname')),'file')==2
            P = load(fullfile(getpref('cellbase','fname')),'PREFERENCES');
            fields = fieldnames(P);
            if ~isempty(fields)
                P=P.PREFERENCES;
                if isempty(varargin) %get all
                    val = P;
                    val.datapath = getpref('cellbase','datapath');
                    val.name = getpref('cellbase','name');
                    val.fname = getpref('cellbase','fname');
                else
                    if isfield(P,varargin{1})
                        val=P.(varargin{1});
                    else %no parameter in PREFERENCES with that name
                        error('getcbpref: no parameter with the specified name found. Try setcbpref() first.');
                    end
                end
            else %no preference file
                import_preferences();
                val = getcbpref(varargin);
            end
            
        else %no cellbase
            error('getcbpref: No Cellbase database found. Check cellbase filename.');
        end
    end
else
    error('getcbpref: Too many input arguments.')
end

