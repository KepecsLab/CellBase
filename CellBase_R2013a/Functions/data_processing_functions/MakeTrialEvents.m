function MakeTrialEvents(cellid, defTE, Events_TTL, Events_TS, varargin)
%% Parse input
p = inputParser;
addRequired(p,'cellid',@(x) ischar(x));
addRequired(p,'defTE',@(x) isstruct(x) || isa(x, 'function_handle'));
addRequired(p,'Events_TTL');
addRequired(p,'Events_TS');
addOptional(p,'TrialStart_State',1,@(x) isnumeric(x) && isscalar(x) && (x > 0));
addOptional(p,'nameTrialStart','TrialStartTimestamp',@(x) ischar(x));
addOptional(p,'rw',true, @(x) islogical(x) || numel(x) == 1);
parse(p,cellid, defTE, Events_TTL, Events_TS, varargin{:});

[subject,session] = cellid2tags(cellid);

%% Check whether to write TrialEvents file
%       do so if either TrialEvents.mat doesn't exist OR rw is true
%       (default)

if ~exist(fullfile(getpref('cellbase','datapath'),subject,session,getpref('cellbase','session_filename')),'file') || p.Results.rw

    tsNeur = Events_TS(Events_TTL==p.Results.TrialStart_State);
    
    if isa(defTE, 'function_handle')
        [TE, tsBhv] = feval(defTE); % expects @defTE to output a struct of 
                                    % InterestingThings (a 1xnTrials vector of event times relative to trial start)
                                    % and tsBhv (trial start times from Bpod)
    elseif isstruct(defTE)  % defTE is a scalar struct of 1xnTrials vectors; ideally should be a 1xnTrials struct array
        tsBhv = p.Results.defTE.(p.Results.nameTrialStart); % Should work for TE = Bpod data file
        TE = rmfield(defTE,p.Results.nameTrialStart);
    end
    
    %% Check whether number of timestamps match between Bpod and recording system
    if numel(tsNeur) ~= numel(tsBhv)
        if numel(tsNeur) > numel(tsBhv)
            iShift = trim(tsNeur,tsBhv);
            tsNeur = tsNeur(iShift:numel(tsBhv)+iShift-1);
        elseif numel(tsNeur) < numel(tsBhv)
            iShift = trim(tsBhv,tsNeur);
            allfields = fieldnames(TE);
            for iField = 1:numel(allfields)
                if numel(TE.(allfields{iField})) == numel(tsBhv)
                    TE.(allfields{iField}) = TE.(allfields{iField})(iShift:numel(tsNeur)+iShift-1);
                end
            end
            tsBhv = tsBhv(iShift:numel(tsNeur)+iShift-1);
        end
    end
    %% Write TE file
    TE.TrialStart = tsNeur;
    
    % Save synchronized 'TrialEvents' file
    save(fullfile(getpref('cellbase','datapath'),subject,session,getpref('cellbase','session_filename')),'-struct','TE')
    % if ~isempty(AlignedRecEvents)
    %     save(fullfile(getpref('cellbase','datapath'),subject,session, 'AlignedRecEvents.mat'),'AlignedRecEvents')
    % end
    
end

%% Finds correct behavior and neuro timestamp alignment when different numbers of trials in each
    function [iMax, delta, rhos] = trim(lon,sho)
        lon = lon(:); sho = sho(:);
        delta = numel(lon)-numel(sho);
        rhos = nan(delta+1,1);
        for i = 1:numel(rhos)
            rhos(i) = corr(diff(sho),diff(lon(i:numel(sho)+i-1)));
        end
        assert(max(rhos)>.99)
        iMax = find(rhos==max(rhos));
    end
end