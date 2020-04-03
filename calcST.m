function STarray = calcST(indHS1,indHS2,fs)

% Find array of times between indHS events. indHS1 is for one limb side,
% indHS2 is for other side. Loop through all indHS1 events and find
% preceding HS2 event

for i = 1:length(indHS1)
    ind = find(indHS2 < indHS1(i),1,'last');
    if isempty(ind)
        STarray(i) = nan;
    else
        STarray(i) = (indHS1(i) - indHS2(ind))/fs;
    end
end