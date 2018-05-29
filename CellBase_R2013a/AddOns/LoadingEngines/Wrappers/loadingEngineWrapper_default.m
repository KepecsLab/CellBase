function [TS, WV] = loadingEngineWrapper_default(cellid)
% Used for loading a cell's spike times and waveforms from the original
% tetrode/shank data when no other function is specified in Cellbase
% Preferences.  
% Checks for 'X.ntt' or 'X.dat' files in session folder, where X is the
% tetrode number.  Uses a default Neuralynx or Trodes loading engine,
% respectively.  Fails otherwise.
% Bypass this script by writing a function that takes cellid
% as an argument and returns the time stamps and waveforms for each spike
% in MClust-readable format. (Identical to functionality of an MClust Loading
% engine, but accepting a cellid instead of a filename.) See LOADTRODEDATA
% description for an example.
% 
% See LOADTRODEDATA(cellid)
% 
% EG 05/2018

%% locate session folder
[animal,session,trode,~] = cellid2tags(cellid);
baseDir = getpref('cellbase', 'datapath');
fs = getpref('cellbase', 'filesep');
sessionDir = [baseDir fs animal fs session];

%% find putative recording files (looks for X.ntt, then X.dat)
putNtt = listfiles(sessionDir, [trode '.ntt']);
putDat = listfiles(sessionDir, [trode '.dat']);

if ~isempty(putNtt)
    warning('Found .ntt files in session folder.  Using default Neuralynx LoadingEngine.');
    data_fname = cellid2fnames(cellid, 'ntt')%;
    [TS, WV] = LoadTT_NeuralynxNT(data_fname);
elseif ~isempty(putDat)
    warning('Found .dat files in session folder.  Using default Trodes LoadingEngine.');
    data_fname = cellid2fnames(cellid, 'SGdat')%;
    [TS, WV] = MClustTrodesLoadingEngine(data_fname);  
else
    error('No loadingEngineWrapper set in Preferences, nor an .ntt or .dat file for this [tet]rode in session folder.')
end

end