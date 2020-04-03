function indHS = getHS(markerPos)

% markerPos is vertical position of marker during period of interest
% (walking forward or backward).

% Get index of vertical markerPos where HS occurred based on extrema. 
x = 1:length(markerPos);
[p,ind] = findpeaks(markerPos,x,'MinPeakProminence',0.01);
indHS = [];
% Find minima between peaks
for i = 2:length(ind)
    [m(i-1),temp] = min(markerPos(ind(i-1):ind(i)));
    indHS(i-1) = temp + ind(i-1) - 1;
end

% For last step, look for min between last peak and end of time
% period of interest
[m(i),temp] = min(markerPos(ind(i):end));
indHS(i) = temp + ind(i) - 1;

% plot(markerPos);hold on,plot(indHS,m,'x'),plot(ind,p,'o');
