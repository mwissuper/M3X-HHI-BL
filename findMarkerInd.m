function indMarker = findMarkerInd(s,subj,markerID)

% Take input string s of name of marker and find which index it is in
% MarkerID. Check for cases of 1 or 2 subjects

% check if two participants
if subj < 10
    str = sprintf('BL0%i:%s',subj,s); % 25 spaces for label name
else
    str = sprintf('BL%i:%s',subj,s);
end
str(end+1:30) = ' '; % add blanks
for i = 1:length(markerID)
    idx(i) = strcmp(markerID(i,:),str);
end
indMarker = find(idx == 1,1,'first');

% check if one participant
if isempty(indMarker)
    str = s; % 25 + 5 spaces for label name
    str(end+1:30) = ' ';
    for i = 1:length(markerID)
        idx(i) = strcmp(markerID(i,:),str);
    end
    indMarker = find(idx == 1,1,'first');
end