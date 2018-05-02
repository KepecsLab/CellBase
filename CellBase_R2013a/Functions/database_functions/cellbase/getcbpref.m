function val = getcbpref(name)
% GETCBPREF() returns current cellbase's preference value for preference name STR.
% Cellbase-specific preferences are stored in cellbase file
% (getpref('cellbase','fname')). Note than cellbase-general preferences
% (cellbase name, data folder, cellbase file) are manipulated with
% getpref('cellbase')/setpref('cellbase',...).
%
% Returns VAL, the value of preference NAME. Returns FALSE if n preference
% with NAME exists.
%
% See also: setcbpref(name,value), default_preferences()
%
% TO 05/2018

if any(strcmp({'datapath','name','fname','filesep','cellbases'},name)) %global settings
    warning('Please use getpref to get global cellbase settings.');
    val = getpref('cellbase',name);
else %correct use for cellbase-specific parameters
    if exist(fullfile(getpref('cellbase','fname')),'file')==2
        P = load(fullfile(getpref('cellbase','fname')),'PREFERENCES');
        fields = fieldnames(P);
        if ~isempty(fields)
            P=P.PREFERENCES;
            if isfield(P,name)
                val=P.(name);
            else %no parameter in PREFERENCES with that name
                val = false;
                warning('getcbpref: no parameter with the specified name found. Try setcbpref() first.');
            end
        else %no preference file
            import_preferences();
            val=getcbpref(name);
        end
        
    else %no cellbase
        error('getcbpref: No Cellbase database found. Check cellbase filename.');
    end
end

end