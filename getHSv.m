function indHS = getHSv(markerPos,B,A,subj,trial,side)

% markerPos is AP position of marker during period of interest
% (walking forward or backward). Get HS event as when AP vel first drops
% below thresh after each peak in v.
% Use filter coeff's calculated previously. % Last flag for L or R for subj
% 6 trial 20 issue of just remove for L foot

vThresh = 0.002;
    
v = filtfilt(B,A,diff(markerPos));
x = 1:length(v);
[p,indP] = findpeaks(v,x,'MinPeakProminence',0.01);
if subj == 6  % special cases, need to find all the steps using custom code
    if trial == 13 && side == 2
        indP = [8 indP];
        p = [v(8); p];
    elseif trial == 20 && side == 1
        p(1) = []; indP(1) = [];
    elseif trial == 27 && side == 2
        indP(end+1) = 968;
        p(end+1) = v(968);
    elseif trial == 30 
        if side == 2
            indP = [53 indP];
            p = [v(53); p];
        end
    elseif trial == 32 && side == 1
        indP(1) = [];
        p(1) = [];
    elseif trial == 33 && side == 2
        indP = [61 indP];
        p = [v(61); p];
    elseif trial == 38 && side == 2
        indP = [82 indP];
        p = [v(82); p];
    elseif trial == 42 && side == 2
        indP = [72 indP];
        p = [v(72); p];
    elseif trial == 45 && side == 2
        indP = [86 indP];
        p = [v(86); p];
    elseif trial == 47 
        if side == 1
            indP(2:3) = [];
            p(2:3) = [];
        else % R side
            indP = [98 indP];
            p = [v(98); p];
        end
    elseif trial == 53 && side == 1
        indP(2) = [];
        p(2) = [];
    end
end
indHS = [];
% Find thresh crossing between peaks
vOffset = [v(2:end); nan];
for i = 2:length(indP)
    if subj == 6 
        if (trial == 47 && side == 2 && i == 2) || (trial == 53 && side == 1 && i == 2) || (trial == 47 && side == 1 && i == 2)
            indHS(i-1) = find(v(indP(i-1):indP(i)) > vThresh & vOffset(indP(i-1):indP(i)) < vThresh,1,'first') + indP(i-1) - 1;
        else
            indHS(i-1) = find(v(indP(i-1):indP(i)) > vThresh & vOffset(indP(i-1):indP(i)) < vThresh,1,'last') + indP(i-1) - 1;
        end
    else
        indHS(i-1) = find(v(indP(i-1):indP(i)) > vThresh & vOffset(indP(i-1):indP(i)) < vThresh,1,'last') + indP(i-1) - 1;
    end
end
if subj == 6 && trial == 43
    indHS(i) = find(v(indP(i):end) > vThresh & vOffset(indP(i):end) < vThresh,1,'first') + indP(i) - 1;
else
    % For last event, look for last thresh crossing before end of trial
    indHS(i) = find(v(indP(i):end) > vThresh & vOffset(indP(i):end) < vThresh,1,'last') + indP(i) - 1;
end

% subplot(2,1,1),plot(markerPos);hold on,plot(indHS,markerPos(indHS),'x'); ylabel('AP pos (m)');
% subplot(2,1,2),plot(v);hold on,plot(indHS,v(indHS),'x'),plot(indP,v(indP));
% ylabel('AP vel (m/s)');



