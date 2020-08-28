function new_cellids = findnewcells()
%FINDNEWCELLS   FIND new cells in CellBase.

% Get CellBase preferences
cellbase_datapath = getpref('cellbase','datapath');

% Find all cells and sort the new ones

all_cellids = findallcells;

old_cellids = listtag('cells');

if isempty(old_cellids)
    new_cellids = all_cellids;  % setdiff fails for empty set
else
    new_cellids = setdiff(all_cellids,old_cellids);
end


