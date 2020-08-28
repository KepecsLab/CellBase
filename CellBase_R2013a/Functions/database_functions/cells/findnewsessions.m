function [new_sessions_rats,new_sessions] = findnewsessions()
%FINDNEWCELLS   FIND new sessions in CellBase (without ANY cell).


newcells = findnewcells();

for f = 1:length(newcells)
    [rat,ses]= cellid2tags(newcells{f});
    Session{f}=ses;
    Rat{f}= rat;
end


[unique_sessions,ii] = unique(Session);
unique_sessions=unique_sessions(:);

unique_session_rats = Rat(ii);

all_sessions = listtag('session');
all_sessions = all_sessions(:,2);

[new_sessions,ii] = setdiff(unique_sessions,all_sessions);

new_sessions_rats = unique_session_rats(ii);