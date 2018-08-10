function setcbpref(name,val)
% SETCBPREF() sets current cellbase's preference value for preference NAME
% to VAL.
% Cellbase-specific preferences are stored in cellbase file
% (getpref('cellbase','fname')). Note that cellbase-general preferences
% (cellbase name, data folder, cellbase file) are manipulated with
% getpref('cellbase')/setpref('cellbase',...).
%
% See also: getcbpref(name), default_preferences()
%
% TO 05/2018, FS 08/2018

if any(strcmp({'datapath','name','fname','filesep','cellbases'},name)) %global preferences
    warning('Please use setpref to set global cellbase settings.');
    setpref('cellbase',name,val);
else %correct use for cellbase-specific preferences
    if exist(fullfile(getpref('cellbase','fname')),'file')==2
        CB = load(fullfile(getpref('cellbase','fname')));
        fields = fieldnames(CB);
        if ~isempty(fields)
            if ischar(val) && (val(1) == '@') % it's to be a function handle % FS MOD
                CB.PREFERENCES.(name)=str2func(val); 
            else
                CB.PREFERENCES.(name)=val;                 
            end
            save(fullfile(getpref('cellbase','fname')),'-struct','CB')
        else %no preference file
            import_preferences();
            setcbpref(name,val);
        end
        
    else %no cellbase
        error('setcbpref: No Cellbase database found. Check cellbase filename.');
    end
end