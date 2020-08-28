function [new_sessions_rats,new_sessions] = findnewcellsessions()
%FINDNEWCELLS   FIND new sessions in CellBase with non-added cells.


newcells = findnewcells();

for f = 1:length(newcells)
    [rat,ses]= cellid2tags(newcells{f});
    Session{f}=ses;
    Rat{f}= rat;
end


[unique_sessions,ii] = unique(Session);
new_sessions=unique_sessions(:);

new_sessions_rats = Rat(ii);
