function [TS, WV] = loadtrodedata(cellid)
% Loads a cell's spike times and waveforms from the original
% tetrode/shank data. Expects a loading engine wrapper written by user and
% specified via setcbpref('LoadingEngineWrapper_fname', 'yourLEWrapper').
% Calls LOADINGENGINEWRAPPER_DEFAULT otherwise.  
% 
% [TS, WV] = yourLEWrapper(cellid) is a function that takes cellid
% as an argument and returns the time stamps and waveforms for each spike
% in the cluster.  The format should match that output by an MClust loading
% engine, hence the name "loadingEngineWrapper".  
%
% EXAMPLE:
% function [TS, WV] = myLoadingEngineWrapper(cellid)
%       trodeDataFN = cellid2fname(cellid, 'ntt');
%       [TS, WV] = myNeuralynxLoadingEngine(trodeDataFN);
% end
% 
% EG 05/2018
    
    loadingEngineWrapper = getcbpref('LoadingEngineWrapper_fname');
    
    [~,prefFN,~] = fileparts(loadingEngineWrapper);
    try strcmp(prefFN, 'loadingEngineWrapper_default')
    catch
        warning('No LoadingEngineWrapper_fname set in Preferences.  Using default.');
    end
 
[TS, WV] = feval(loadingEngineWrapper, cellid);   
end

    