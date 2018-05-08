function import_preferences()
% IMPORT_PREFERENCES creates a struct PREFERENCES containing default parameters
% or loads preferences from getpref('cellbase'). This is a
% convenience function to maintain backwards compatibility to cellbase
% 5fbf007 (12 Sep 2017) and prior.

% See also default_preferences(), setcbpref(name,value), getcbpref(name), initcb()
%
% TO 05/2018

str = 'No PREFERENCE variable found. Likely due to use of older cellbase version. Please choose if you want to use default parameters or import from getpref("cellbase"). ';
str=strcat(str,'Preferences can be manipulated using getcbpref and setcbpref and are stored in your cellbase file.');

A = questdlg(str,'convert_preferences','Default','Import','Default');

PREFERENCES = default_preferences();

switch A
    case 'Default'
        %
    case 'Import'
        PREFERENCES.TrialEvents_fname = getpref('cellbase','TrialEvents_filename');
        PREFERENCES.StimEvents_fname = getpref('cellbase','StimEvents_filename');
        PREFERENCES.Spikes_cell_pattern = getpref('cellbase','cell_pattern');
        PREFERENCES.Spikes_timefactor = getpref('cellbase','timefactor');
end

CB = load(fullfile(getpref('cellbase','fname')));
CB.PREFERENCES=PREFERENCES;
save(fullfile(getpref('cellbase','fname')),'-struct','CB')
