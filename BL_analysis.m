% Plot time series data for an individual subject. Calculate metrics for
% participant. Do all filtering in Nexus, so just look at
% rawData.video.markers for marker data

clear; clc; close all;

% Analysis options

subj = 9; % Total subj's 1:11. Do not trust FT1 for BL04 (looks not connected). BL01 force data might also be bad but need to check
analyzeForces = 1; % Tension is positive. Not confident of directions of off-axis force vectors unless use recoverForces so don't try to interpret directionality
fcA = 10; % cutoff freq for force in hz
fcM = 10; % cutoff freq for marker data in hz

%% Plotting options
plotTorso = 0; % 1 plot AP dir to check speed calculations
plotFootHS = 0; % Check HS events are ID'ed correctly
plotVyTorso = 0; % Check begninning of moving vs. not moving
plotCheckSLST = 0;
plotFzV = 0; % Plot raw voltage for Fz for all trials as check on bias
plotFz = 0; % Plot force in N with HS events
% Rarely used
plotArmLen = 0; % Plot arm lengths of each participant in AP dir vs. time
plotVyFz = 0; % Check what happened with force in plots with strange velocity (non-constant) during middle 4 steps
plotCorr = 0; % Plot correlation of FT1 and FT2 Fz
plotFxy = 0;
plotFzBias = 0; % Plot Fz volt vs. time for each axis for bias voltage trial

if plotVyTorso == 1 || plotFzV == 1 || plotTorso == 1 || plotFootHS == 1
    plotAllTrials = 1;
else
    plotAllTrials = 0; % plot only partner trials
end

pcount = 0;

%% Constants

if subj == 2
    trialArray = 1:56; %[1:49 51 53:56]; % Trial 50 and 52 participant walked out of capture volume during last few forward steps, but if analysis is of middle 4 steps, then ok to use. Trial 50 won't export, so skip it.
elseif subj == 3
    trialArray = [1:8 10:21 23:26 28 29 31:56 58 59]; % Some bad trials (22, 27, 30), trial 9 started with left foot and didn't catch it during expt. Incorrectly repeated pref in trial 57 instead of asymm
elseif subj == 5
    trialArray = 1:59; % Trial 60 didn't record for some reason
elseif subj == 6
    trialArray = [1:11 13:61]; % Trial 12 poor fill feet
elseif subj == 9
    trialArray = [1:25 27:61]; % Skip bad trial 26
elseif subj == 10
    trialArray = [1:18 20:61]; % Skip bad trial
elseif subj == 11
    trialArray = [1:27 29:60]; % Skip bad trial 28 kicked MH forward walk. Trial 61 didn't record!
else
    trialArray = 1:60;
end
sourcefolder = cd;
if subj < 10
    subjFolder = sprintf('BL0%i',subj);
    configFile = sprintf('BL0%i_worksheet.mat',subj);
else
    subjFolder = sprintf('BL%i',subj);
    configFile = sprintf('BL%i_worksheet.mat',subj);
end
 
temp = load(configFile);
config = temp.data;

colors(1,:) = [0.00,0.45,0.74]; % nice blue
colors(2,:) = [0.85,0.33,0.10]; % nice red
colors(3,:) = [0.47,0.67,0.19]; % nice green

markerArray = {'CLAV','C7','LHEE','RHEE'}; % Markers to be analyzed

badTrial = zeros(size(trialArray)); % Init to no bad trials

%% First get voltage bias for zero for FT1 for each participant. Use first solo walking trial after or last solo walking trial before partner trials 
if analyzeForces == 1
    if subj == 2
        biasTrial = 45;
    elseif subj == 3
        biasTrial = 45; 
    elseif subj == 8
        biasTrial = 8; % FT1 looks disconnected after trial 48
    else % All subj's after BL04 have same number of trials. Subj 4 not usable due to sensor error throughout.
        biasTrial = 49;
    end
    
    if biasTrial < 10
        filename = sprintf('Trial0%i.mat',biasTrial);
    else
        filename = sprintf('Trial%i.mat',biasTrial);
    end
    trialData = load([sourcefolder,'\',subjFolder,'\',filename]);
    s = size(trialData.rawData.analog.otherid);
    fsM = trialData.VideoFrameRate;  
    % Set up filter for force data
    fsA = trialData.AnalogFrameRate;
    wn = fcA/(fsA/2);
    [Bf,Af] = butter(4,wn,'low');
    
    % Loop through all id names for analog input channels to pull out forces 
    for row = 1:s(1)
        if strcmp(trialData.rawData.analog.otherid(row,:),'FT1 Electric Potential.Fx       ')
            chan.FT1.Fx = row;
        elseif strcmp(trialData.rawData.analog.otherid(row,:),'FT1 Electric Potential.Fy       ')
            chan.FT1.Fy = row;
        elseif strcmp(trialData.rawData.analog.otherid(row,:),'FT1 Electric Potential.Fz       ')
            chan.FT1.Fz = row;
        end
    end
    
    % Filter the force components
    FT1.Fx = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fx));
    FT1.Fy = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fy));
    FT1.Fz = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fz));
    
    % change sampling to match marker data
    Fx1 = downsample(FT1.Fx,fsA/fsM); 
    Fy1 = downsample(FT1.Fy,fsA/fsM);
    Fz1 = downsample(FT1.Fz,fsA/fsM); 

    indStart = 1; indEnd = 120; % Look at short time period before start walking and keep consistent across all trials

    % Just use mean during period of interest as metric
    biasFx1 = nanmean(Fx1(indStart:indEnd)); % voltage
    biasFy1 = nanmean(Fy1(indStart:indEnd));
    biasFz1 = nanmean(Fz1(indStart:indEnd));
%     biasFz1 = -0.1164; % Mean of all individ subj zeros
    % Plot each force signal vs. time to make sure mean is taken during period
    % of relatively little motion

    if plotFzBias == 1
        plotind = 0;

        plotind = plotind + 1;
        subplot(3,1,plotind);
        plot(trialData.mtime(indStart:indEnd),Fz1(indStart:indEnd),'g');
        hline(biasFz1,'g');
        titlename = sprintf('BL%i Trial %i %s',subj,biasTrial,config{biasTrial,2});title(titlename);
        ylabel('Fz1 (V)');

        plotind = plotind + 1;
        subplot(3,1,plotind);
        plot(trialData.mtime(indStart:indEnd),Fx1(indStart:indEnd),'g');
        hline(biasFx1,'g');
        ylabel('Fx1 (V)');

        plotind = plotind + 1;
        subplot(3,1,plotind);
        plot(trialData.mtime(indStart:indEnd),Fy1(indStart:indEnd),'g');
        hline(biasFy1,'g');
        ylabel('Fy1 (V)');
    end
end

%% algorithm parameters for main analysis
vyThresh = 0.1; % vel thresh to start trial (m/s)

if subj >= 4 % now have baseline (eyes open) trials
    numrows = 6;
    numcols = 9;
else
    if plotAllTrials == 1 % plot all trials
        numrows = 8;
        numcols = 6;
    else % plot only experimental trials with partner
        numrows = 8;
        numcols = 4;
    end
end
plotind = 0;
plotind1 = 0; plotind2 = 0; plotind3 = 0; plotind4 = 0;

% Initialize empty arrays so can add to them in loop. Afterwards,
% concatenate them and take means
speedDec = []; speedInc = []; speedPref = []; speedFollow = [];
speedSlowPost = []; speedFastPost = []; speedPrefPost = [];
speedPrefPre = []; speedBase = [];

LSLDec = []; LSLInc = []; LSLPref = []; LSLFollow = [];
LSLSlowPost = []; LSLFastPost = []; LSLPrefPost = [];
LSLPrefPre = []; LSLBase = [];
RSLDec = []; RSLInc = []; RSLPref = []; RSLFollow = [];
RSLSlowPost = []; RSLFastPost = []; RSLPrefPost = [];
RSLPrefPre = []; RSLBase = [];
SLDec = []; SLInc = []; SLPref = []; SLFollow = [];
SLSlowPost = []; SLFastPost = []; SLPrefPost = [];
SLPrefPre = []; SLBase = [];

LcadDec = []; LcadInc = []; LcadPref = []; LcadFollow = [];
LcadSlowPost = []; LcadFastPost = []; LcadPrefPost = [];
LcadPrefPre = []; LcadBase = [];
RcadDec = []; RcadInc = []; RcadPref = []; RcadFollow = [];
RcadSlowPost = []; RcadFastPost = []; RcadPrefPost = [];
RcadPrefPre = []; RcadBase = [];
cadDec = []; cadInc = []; cadPref = []; cadFollow = [];
cadSlowPost = []; cadFastPost = []; cadPrefPost = [];
cadPrefPre = []; cadBase = [];

meanFz1Dec = []; meanFz1Inc = []; meanFz1Pref = []; meanFz1Follow = [];
meanFz1SlowPost = []; meanFz1FastPost = []; meanFz1PrefPost = [];
meanFz1PrefPre = []; meanFz1Base = [];
meanFz2Dec = []; meanFz2Inc = []; meanFz2Pref = []; meanFz2Follow = [];
meanFz2SlowPost = []; meanFz2FastPost = []; meanFz2PrefPost = [];
meanFz2PrefPre = []; meanFz2Base = [];

corrFzDec = []; corrFzInc = []; corrFzPref = []; corrFzFollow = [];
corrFzSlowPost = []; corrFzFastPost = []; corrFzPrefPost = [];
corrFzPrefPre = []; corrFzBase = [];

meanFxy1Dec = []; meanFxy1Inc = []; meanFxy1Pref = []; meanFxy1Follow = [];
meanFxy1SlowPost = []; meanFxy1FastPost = []; meanFxy1PrefPost = [];
meanFxy1PrefPre = []; meanFxy1Base = [];
meanFxy2Dec = []; meanFxy2Inc = []; meanFxy2Pref = []; meanFxy2Follow = [];
meanFxy2SlowPost = []; meanFxy2FastPost = []; meanFxy2PrefPost = [];
meanFxy2PrefPre = []; meanFxy2Base = [];

%% Loop through all trials for participant

for n = 1:length(trialArray)
    clear Markers
    trial = trialArray(n)
    if trial < 10
        filename = sprintf('Trial0%i.mat',trial);
    else
        filename = sprintf('Trial%i.mat',trial);
    end
    trialData = load([sourcefolder,'\',subjFolder,'\',filename]);
    if isfield(trialData,'VideoFrameRate')
        fsM = trialData.VideoFrameRate; 
    else
        fsM = trialData.data.video.samplerate;
    end 
    % Set up filter for marker data
    wn = fcM/(fsM/2);
    [Bm,Am] = butter(4,wn,'low');
    Markers = trialData.rawData.video.markers;
    
    % c3d files put in zeros for
    % gaps so need to find these and make these nan's!
    % Loop through all markers
    % Will only work for markers that are mostly labeled, so choose which
    % ones to preprocess this way.
    for i = 1:length(markerArray)
        temp = findMarkerInd(markerArray{i},subj,trialData.MarkerID);
        if ~isempty(temp) % marker not present, sometimes it's the case for C7
            indMarker(i) = temp;
            msg = sprintf('%s not present',markerArray{i});
        else
            indMarker(i) = nan;
        end
    end
    clear temp
        
    for i = indMarker
        if ~isnan(i)
            Markers(:,i,:) = replaceZeros(trialData.rawData.video.markers(:,i,:));
        end
%         if isempty(find(isnan(Markers(:,i,:)),1,'first')) % no nan values, ok to use filtfilt
%             Markers(:,i,:) = filtfilt(Bm,Am,Markers(:,i,:)); % will error out if there are nan's!
%         else
%             Markers(:,i,:) = filtIgnoreNan(Bm,Am,order,Markers(:,i,:));
%         end
    end

    cond = config{trial,2};
    if ~strcmp(cond,'Asymmetric')
        
        %% Find when forward walking started and ended based on clav marker
        
        % Find index of MarkerID for particular marker
        temp.indClav = findMarkerInd('CLAV',subj,trialData.rawData.video.markersid);
        clav = squeeze(Markers(:,temp.indClav,:));
        temp = [];
        clav = clav./1000; % convert to m
        
        % If gaps in clav, that means significant gap and can't use this
        % marker. Use C7 in this case.
        if ~isempty(find(isnan(clav),1,'first')) % sig. gaps
            temp.indC7 = findMarkerInd('C7',subj,trialData.rawData.video.markersid);
            c7 = squeeze(Markers(:,temp.indC7,:))./1000;
            temp = [];
            torso = c7;
            disp('c7')
        else
            torso = clav;
        end
        
        if subj == 6 && trial == 53 % special case of gap is at end of trial beyond steps we care about
            torso = clav;
            disp('clav');
        end
        
        vyTorso = diff(torso(:,2))*fsM; % (m/s)
        vyTorsoOffset = [vyTorso(2:end); nan];
        indStart = find(vyTorso > vyThresh,1,'first');
        indEnd = find(vyTorso(indStart:end) > vyThresh & vyTorsoOffset(indStart:end) < vyThresh,1,'last') + indStart;

        %% Plot torso marker pos to check that beginning and end of fwd walking period found correctly
        if plotTorso == 1
            plotind = plotind + 1;
            subplot(numcols,numrows,plotind)
%             yyaxis left
            plot(trialData.mtime,torso(:,2)),ylabel('Torso AP pos (m)')
%             yyaxis right
%             plot(trialData.mtime(2:end),vyTorso),ylabel('Torso AP vel (mm/s)'); hold on;
            vline([trialData.mtime(indStart) trialData.mtime(indEnd)],'k-')
            xlim([trialData.mtime(indStart)-0.2 trialData.mtime(indEnd)+0.2])
            if plotind == 1
                titlename = sprintf('BL%i Trial %i',subj,trial); 
            else
                titlename = sprintf('Trial %i',trial); 
            end
            title(titlename);
        end
        
        %% Find HS events from heel marker data. Only use ones from forward walking period.
        % Use ankle marker for BL06 since heel coverage bad.
        % Find index of MarkerID for particular marker
       
        if subj == 6
            temp.indLfoot = findMarkerInd('LANK',subj,trialData.rawData.video.markersid);
            temp.indRfoot = findMarkerInd('RANK',subj,trialData.rawData.video.markersid);
        else
            temp.indLfoot = findMarkerInd('LHEE',subj,trialData.rawData.video.markersid);
            temp.indRfoot = findMarkerInd('RHEE',subj,trialData.rawData.video.markersid);
        end
        Lfoot = squeeze(Markers(:,temp.indLfoot,:));
        Lfoot = Lfoot./1000; % convert to m
        Rfoot = squeeze(Markers(:,temp.indRfoot,:));
        Rfoot = Rfoot./1000; % convert to m

        if subj == 6
            indLHS = getHSv(Lfoot(1:indEnd,2),Bm,Am,subj,trial,1);
            indRHS = getHSv(Rfoot(1:indEnd,2),Bm,Am,subj,trial,2);
        else
            indLHS = getHS(Lfoot(1:indEnd,3));
            indRHS = getHS(Rfoot(1:indEnd,3));
        end
        
        if ~(subj==6 && (trial == 3 || trial == 43 || trial == 55 || trial == 56)) % special case where need this LHS (for SLST calc's) and vyTorso is poor due to marker issues
            indLHS(find(indLHS <= indStart | indLHS >= indEnd)) = [];
        end
        
%         indRHS(find(indRHS <= indStart | indRHS >= indEnd)) = [];

        temp = [];
                
        %% Exceptions where algo didn't work and manually correct event 
        % based on visual observation in Nexus. Way this is written with
        % elseif, must remove all events for a trial at once
        if subj == 2 
            if trial == 14
                indRHS(3) = 624;
            elseif trial == 18
                indRHS(2) = 356;
                indRHS(3) = 559;
            elseif trial == 27
                indRHS(4) = 1120;
            elseif trial == 32
                indLHS(4) = 1066;
            elseif trial == 38
                indRHS(2) = [];
            end
        elseif subj == 3
            if trial == 5
                indRHS(2) = 408;
                indRHS(3) = 648;
                %indLHS(2) = 281;
            elseif trial == 15
                indRHS(2) = 381;
            end
        elseif subj == 4
            if trial == 4
                indRHS(4) = 624;
                indRHS(5) = 764;
            elseif trial == 5
                indRHS(4) = 756;
                indLHS(end) = 846;
            elseif trial == 9
                indRHS(2) = 514;
            elseif trial == 14
                indLHS(end) = 710;
            elseif trial == 16
                indLHS(end) = 764;
            elseif trial == 24
                indRHS(3) = 743;
                indLHS(3) = 865;
            elseif trial == 27
                indLHS(end-1) = [];
            elseif trial == 31
                indLHS(3) = [];
                indRHS(2) = 425;
                indRHS(4) = 1082;
            elseif trial == 34
                indLHS(3) = [];
            elseif trial == 39
                indLHS(end) = 732;
            elseif trial == 47
                indRHS(3) = 595;
            end
        elseif subj == 5
            if trial == 1
                indRHS = [164 indRHS];
            elseif trial == 3
                indRHS = [176 indRHS];
            elseif trial == 22
                indRHS = [198 indRHS];
            elseif trial == 40
                indRHS(3) = [];
            end
        elseif subj == 6
            if trial == 17
                indRHS = [133 indRHS(1:2) 720 indRHS(3)];
            elseif trial == 18
                indLHS = indLHS(end-3:end);
            elseif trial == 21
                indLHS(1:2) = [];
            elseif trial == 22 || trial == 24
                indLHS(1:6) = [];
            elseif trial == 25
                indRHS([1 3]) = [];
                indLHS(1) = [];
            elseif trial == 27
                indLHS(1:2) = [];
            elseif trial == 31
                indLHS(2:3) = [];
            elseif trial == 34
                indRHS(1:3) = [];
            elseif trial == 36 || trial == 45
                indLHS(1:6) = [];
            elseif trial == 41
                indRHS(1:2) = [];
                indLHS(1:5) = [];
            elseif trial == 46
                indRHS(1:2) = [];
            elseif trial == 57
                indLHS(1:2) = [];
            end
        elseif subj == 7
            if trial == 14
                indRHS(3) = 1372-629;
            elseif trial == 21
                indRHS(2:3) = [846 1091]-400;
                indLHS(1:3) = [725 958 1212] - 400;
            elseif trial == 39
                indRHS(3:4) = [693 943];
            end
        elseif subj == 8
            if trial == 2
                indRHS(4) = [];
            elseif trial == 7
                indLHS(3) = [];
            end
        elseif subj == 9
            if trial == 3
                indRHS([1 3 6]) = [];
            elseif trial == 6 || trial == 8 || trial == 14
                indRHS(4) = [];
            elseif trial == 7
                indRHS([1 4 6]) = [];
                indLHS(4) = [];
            elseif ~isempty(find(trial == [11:13 17 31 33 50 55:56],1,'first'))
                indRHS(1) = [];
            elseif trial == 36 
                indLHS(4) = [];
            elseif trial == 53
                indRHS(3) = [];
            elseif trial == 54
                indRHS([1 4]) = [];
            elseif trial == 60
                indLHS(4) = [];
                indRHS(3) = [];
            end
        elseif subj == 11
            if trial == 27
                indRHS(4) = 1087;
            end
        end
        
        %% Update start and end period of analysis to HS events of interest
        indStart = indRHS(2); indEnd = indRHS(4);
                
        % Check that torso vy is constant during this period by fitting a
        % line to vy during this period and see if regression is
        % statistically significant. Found many trials have statistically
        % significant slope, so choose a threshold for slope.
%         x = [(indStart:indEnd)' ones(size(vyTorso(indStart:indEnd)))];
        clear mdl p m
        mdl = fitlm(indStart:indEnd,vyTorso(indStart:indEnd));
        p = coefTest(mdl);
        if p < 0.05 % If linear regression is significant, check the slope. If slope is larger than threshold, mark the trial to throw out from final analysis
            temp = mdl.Coefficients;
            tempm = temp.Estimate(:,1); % get coeff's for plotting later. 
            m(1) = tempm(2); m(2) = tempm(1); % need to reverse order for polyval. m(1) is slope
            if abs(m(1)) > 0.001 % m/s^2
                badTrial(n) = 1;
            end
            clear temp tempm
        end
        
        %% Plot vel of torso marker to check that beginning and end of fwd 
        % walking period found correctly. Plot HS events of period to see
        % if changing period of interest to middle HS events results in
        % picking out constant vel period.
        if plotVyTorso == 1
            plotind = plotind + 1;
            subplot(numcols,numrows,plotind)
            plot(trialData.mtime(2:end),vyTorso),ylabel('Torso AP vel (m/s)'); hold on;
            % xlim([trialData.mtime(indStart)-0.2 trialData.mtime(indEnd)+0.2])
            ylim([0 1.5])
%             vline([trialData.mtime(indLHS(1)) trialData.mtime(indLHS(end))],'b-')
            % hline(vyThresh,'k--');
            % If bad trial (sig linear slope), plot slope
            if badTrial(n) == 1
                plot(trialData.mtime(indStart:indEnd),polyval(m,indStart:indEnd),'r--');
            end
            if plotind == 1
                titlename = sprintf('BL%i Trial %i %s',subj,trial,cond); 
            else
                titlename = sprintf('Trial %i %s',trial,cond); 
            end
            title(titlename); vline([trialData.mtime(indStart) trialData.mtime(indEnd)],'k-')
        end
        %% Calculate avg speed from distance traveled by torso marker during
        % forward walking portion
        dist(n) = torso(indEnd,2) - torso(indStart,2);
        speed(n) = dist(n)/(trialData.mtime(indEnd)-trialData.mtime(indStart));
        
        %% Plot marker pos to check HS found correctly. Threw away first 
        % init. step and last collection step, so should see 7 steps.
        if plotFootHS == 1
            numcols = 4; numrows = 3;
            plotind = plotind + 1;
            subplot(numcols,numrows,plotind),
            hold on;
            plot(trialData.mtime,Lfoot(:,3),'color',colors(1,:));
            plot(trialData.mtime,Rfoot(:,3),'color',colors(2,:));
            plot(trialData.mtime(indLHS),Lfoot(indLHS,3),'x','color',colors(1,:));
            plot(trialData.mtime(indRHS),Rfoot(indRHS,3),'x','color',colors(2,:));
            ylabel('Foot vert pos (m)')
            if plotind == 1
                h = legend('Left','Right','LHS','RHS');
                set(h,'orientation','horizontal');
            end
            vline([trialData.mtime(indStart) trialData.mtime(indEnd)],'k-')
            % xlim([trialData.mtime(indStart)-0.2 trialData.mtime(indEnd)+0.2])
            titlename = sprintf('Trial %i',trial); title(titlename);
        end
        
        %% Calculate SL for L and R separately and combined. Only use values within period of interest
        
        LSLarray = calcSL(Lfoot(:,2),Rfoot(:,2),indLHS(2:3));
        RSLarray = calcSL(Rfoot(:,2),Lfoot(:,2),indRHS(2:4));
        SLarray = [LSLarray; RSLarray];
        LSLtrial(n) = nanmean(LSLarray);
        RSL(n) = nanmean(RSLarray);
        SL(n) = nanmean(SLarray);
                
        %% Calculate ST and cadence for L and R separately and combined
        % Only use period of interest
        
        LSTarray = calcST(indLHS(2:3),indRHS(2:4),fsM);
        RSTarray = calcST(indRHS(2:4),indLHS(2:3),fsM);
        STarray = [LSTarray RSTarray];
        LST(n) = nanmean(LSTarray);
        RST(n) = nanmean(RSTarray);
        ST(n) = nanmean(STarray);
        
        Lcad(n) = 1./LST(n).*60; % convert to steps per min.
        Rcad(n) = 1./RST(n).*60; 
        cad(n) = 1./ST(n).*60;
        
        %% Plot to check SL and ST calculated correctly
        % Plot AP position L and R feet (heel or ankle)
        % Plot SL at LHS and RHS events of interest using LSL and RSL
        % Plot timing of LHS and RHS using LST and RST
        if plotCheckSLST == 1
            numcols = 6; numrows = 4;
            plotind = plotind + 1;
            subplot(numrows,numcols,plotind);
            plot(trialData.mtime,Lfoot(:,2),'b'),hold on;
            plot(trialData.mtime,Rfoot(:,2),'r'),hold on;
            if plotind == 1
                legend('L','R','LHS','RHS');
            end
            ind1 = indRHS(2:3); ind2 = indLHS(1:3);
            % Remove nan terms for plotting
            x = trialData.mtime(ind1)+LSTarray;
            y = ind1+time2ind(LSTarray,fsM);
            y(isnan(x)) = []; LSLarray(isnan(x)) = []; x(isnan(x)) = []; 
            x(isnan(y)) = []; LSLarray(isnan(y)) = []; y(isnan(y)) = []; % remove nan's from ind1 or LSTarray
            x(isnan(LSLarray)) = []; y(isnan(LSLarray)) = []; LSLarray(isnan(LSLarray)) = []; % remove nan's from LSLarray
            % Plot Lfoot
            plot(x,Rfoot(y,2)+LSLarray,'bx');
            % Remove nan's
            x = trialData.mtime(ind2)+RSTarray;
            y = ind2+time2ind(RSTarray,fsM);
            y(isnan(x)) = []; RSLarray(isnan(x)) = []; x(isnan(x)) = []; 
            x(isnan(y)) = []; RSLarray(isnan(y)) = []; y(isnan(y)) = []; 
            x(isnan(RSLarray)) = []; y(isnan(RSLarray)) = []; RSLarray(isnan(RSLarray)) = []; 
            % Plot Rfoot
            plot(x,Lfoot(y,2)+RSLarray,'rx');
            box off,set(gca,'tickdir','out');
            xlabel('Time (s)');ylabel('AP Heel position');
            vline([trialData.mtime(indStart) trialData.mtime(indEnd)],'k-')
            titlename = sprintf('Trial %i',trial); title(titlename);
        end
        
        %% Force calculations. 
        if analyzeForces == 1 && badTrial(n) ~= 1
            if plotFzV == 1 % Plot raw voltages (no filtering) for each trial. A way to check that bias voltage value makes sense.
                s = size(trialData.rawData.analog.otherid);
                % Loop through all id names for analog input channels to pull out forces 
                for row = 1:s(1)
                    if strcmp(trialData.rawData.analog.otherid(row,:),'FT1 Electric Potential.Fx       ')
                        chan.FT1.Fx = row;
                    elseif strcmp(trialData.rawData.analog.otherid(row,:),'FT1 Electric Potential.Fy       ')
                        chan.FT1.Fy = row;
                    elseif strcmp(trialData.rawData.analog.otherid(row,:),'FT1 Electric Potential.Fz       ')
                        chan.FT1.Fz = row;
                    end
                end
                % Filter force data (this is the only time it's filtered)
                FT1.Fx = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fx));
                FT1.Fy = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fy));
                FT1.Fz = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fz));
                
                % change sampling to match marker data
                Fx1 = downsample(FT1.Fx,fsA/fsM); 
                Fy1 = downsample(FT1.Fy,fsA/fsM);
                Fz1 = downsample(FT1.Fz,fsA/fsM); 

                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                plot(trialData.mtime,Fz1,'g');
                hline(biasFz1,'g'); vline([trialData.mtime(indStart) trialData.mtime(indEnd)],'k-');
                if trial == biasTrial % Check that biasVoltage matches mean during this period
                    vline(trialData.mtime(120),'r-'); % Time up to this point used for fzBias calc
                end
                if plotind == 1
                    titlename = sprintf('BL%i T%i %s',subj,trial,config{trial,2});
                else
                    titlename = sprintf('T%i %s',trial,config{trial,2});
                end
                if subj == 9
                    ylim([-0.2 0.05]);
                elseif subj == 11
                    ylim([-0.15 0.055]);
                end
%                 ylim([-0.45 0.05]); 
                box off; set(gca,'tickdir','out');
                title(titlename);
                ylabel('Fz1 (V)'); xlabel('Time (s)');
            end
            if ~strcmp(cond,'Baseline') & ~strcmp(cond,'PreferredPre') & ~strcmp(cond,'SlowPost') & ~strcmp(cond,'FastPost') & ~(subj == 5 && trial == 22)% weird trial no force data
                 s = size(trialData.rawData.analog.otherid);
                 % Loop through all id names for analog input channels to pull out forces 
                 for row = 1:s(1)
                     if strcmp(trialData.rawData.analog.otherid(row,:),'FT1 Electric Potential.Fx       ')
                         chan.FT1.Fx = row;
                     elseif strcmp(trialData.rawData.analog.otherid(row,:),'FT1 Electric Potential.Fy       ')
                         chan.FT1.Fy = row;
                     elseif strcmp(trialData.rawData.analog.otherid(row,:),'FT1 Electric Potential.Fz       ')
                         chan.FT1.Fz = row;
                     elseif strcmp(trialData.rawData.analog.otherid(row,:),'FT2 Electric Potential.Fx       ')
                         chan.FT2.Fx = row;
                     elseif strcmp(trialData.rawData.analog.otherid(row,:),'FT2 Electric Potential.Fy       ')
                         chan.FT2.Fy = row;
                     elseif strcmp(trialData.rawData.analog.otherid(row,:),'FT2 Electric Potential.Fz       ')
                         chan.FT2.Fz = row;
                     end
                 end
                 
                 % Filter force data (this is the only time it's filtered)
                 % and subtract out voltage bias when known
                 FT1.Fx = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fx))-biasFx1;
                 FT1.Fy = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fy))-biasFy1;
                 FT1.Fz = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fz))-biasFz1;
                 FT2.Fx = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT2.Fx));
                 FT2.Fy = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT2.Fy));
                 FT2.Fz = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT2.Fz));
                
                 % change sampling to match marker data
                 Fx1 = downsample(FT1.Fx,fsA/fsM); 
                 Fy1 = downsample(FT1.Fy,fsA/fsM);
                 Fz1 = downsample(FT1.Fz,fsA/fsM); 
                 Fx2 = downsample(FT2.Fx,fsA/fsM); 
                 Fy2 = downsample(FT2.Fy,fsA/fsM);
                 Fz2 = downsample(FT2.Fz,fsA/fsM);

                 % Lump together off-axis forces bc unsure of marker
                 % orientation without recoverForces code and bc markers
                 % possible incorrectly attached during BL01-BL03 or BL04
                 Fxy1 = sqrt(Fx1.^2 + Fy1.^2);
                 Fxy2 = sqrt(Fx2.^2 + Fy2.^2);

                 % Calculate correlation FT1 and FT2 voltages after
                 % filtering
                 [r,pval(n)] = corr(Fz1(indStart:indEnd),Fz2(indStart:indEnd));
                 if pval(n) < 0.05 % if corr is sig
                     rho(n) = r;
                 else
                     rho(n) = nan;
                 end
                 
                 % Convert to N after filtering. Multiply by factor to
                 % convert to N from V. FT2 gains changed 12/5/19 from 5V
                 % FS to 10V FS.
                 if subj < 10
                     SFz2 = 1/0.01;
                     SFxy2 = 1/0.04;
                 else
                     SFz2 = 1/0.02;
                     SFxy2 = 1/0.08;
                 end
                 SFz1 = 1/0.02;
                 SFxy1 = 1/0.08;
                 Fz1 = Fz1*SFz1; 
                 Fz2 = Fz2*SFz2;
                 Fxy1 = Fxy1*SFxy1;
                 Fxy2 = Fxy2*SFxy2;
                 
                 % Use mean during period of interest as metric
                 meanFz1(n) = nanmean(Fz1(indStart:indEnd));
                 meanFxy1(n) = nanmean(Fxy1(indStart:indEnd));
                 meanFz2(n) = nanmean(Fz2(indStart:indEnd));
                 meanFxy2(n) = nanmean(Fxy2(indStart:indEnd));
                 
                 %% Plot each condition in separate plot in chrono order
                 if ~strcmp(cond,'PreferredPost') && (plotFz == 1 || plotFxy == 1)
                    numrows = 4; numcols = 2;
                    % Choose which figure to plot to depending on cond
                    if strcmp(cond,'Increase')
                        figure(1); plotind1 = plotind1 + 1;
                        plotind = plotind1;
                    elseif strcmp(cond,'Decrease')
                        figure(2); plotind2 = plotind2 + 1;
                        plotind = plotind2;
                    elseif strcmp(cond,'Follow')
                        figure(3); plotind3 = plotind3 + 1;
                        plotind = plotind3;
                    elseif strcmp(cond,'Preferred')
                        figure(4); plotind4 = plotind4 + 1;
                        plotind = plotind4;
                    end
                    subplot(numrows,numcols,plotind),
                    hold on;
                    if plotFz == 1 
                        if plotCorr == 1 % subtract out means in order to look at correlation of 2 signals
                            plot(trialData.mtime,Fz1-mean(Fz1),'color',colors(1,:)),hold on
                            plot(trialData.mtime,Fz2-mean(Fz2),'color',colors(2,:)); % Participant L hand FT2, R hand FT1
                        else
                            plot(trialData.mtime,Fz1,'color',colors(1,:)),hold on
                            plot(trialData.mtime(indStart:indEnd),meanFz1(n)*ones(size(indStart:indEnd)),'--','color',colors(1,:));
%                             plot(trialData.mtime,Fz2,'color',colors(2,:)); % Participant L hand FT2, R hand FT1
                        end
                        if plotind == 1 && strcmp(cond,'Decrease')
                            legend('FT1','FT1 mean');
                        end
                        ylabel('F axial (N)');
                        if subj == 2
                            ylims = [-12 8];
                        elseif subj == 3
                            ylims = [-13 6.5];
                        elseif subj == 4
                            ylims = [-6.5 6];
                        elseif subj == 5
                            ylims = [-20 10];
                        elseif subj == 7
                            ylims = [-14 10];
                        elseif subj == 8
                            ylims = [-4.5 5];
                        elseif subj == 9
                            ylims = [-4.5 4];
                        elseif subj == 10
                            ylims = [-3 9];
                        elseif subj == 11
                            ylims = [-3 8];
                        end
                    elseif plotFxy == 1
                        plot(trialData.mtime,Fxy1);%,trialData.mtime,Fxy2); % Participant L hand FT2, R hand FT1
                        ylabel('F off-axis (N)');
                        if subj == 3
                            ylims = [-30 15];
                        elseif subj == 4
                            ylims = [-4 6];
                        end
                    end
                    ylim(ylims);
                    hline(0,'k-');
                    xlabel('Time (s)')
                    if plotind == 1 
                        titlename = sprintf('BL%i, Tr. %i %s, rho = %.2f',subj,trial,cond,rho(n)); 
                    else
                        titlename = sprintf('Tr. %i, rho = %.2f',trial,rho(n)); 
                    end
                    title(titlename);
                    v = vline(trialData.mtime(indLHS),'-'); set(v,'color',colors(1,:)); 
                    v = vline(trialData.mtime(indRHS),'-'); set(v,'color',colors(2,:)); % Participant L foot red, R foot blue
                    vline([trialData.mtime(indStart) trialData.mtime(indEnd)],'k-')
                 end
                 
                 %% Plot to check Fz in trials where vy torso is strange
                 if plotVyFz == 1
                     yyaxis left
                     plotind = plotind + 1;
                     subplot(2,2,plotind);
                     plot(trialData.mtime,Fz1,'color',colors(1,:)),hold on
                     plot(trialData.mtime(indStart:indEnd),meanFz1(n)*ones(size(indStart:indEnd)),'--','color',colors(1,:));
                     set(gca,'ycolor',colors(1,:));
                     ylabel('Fz (N)');
                     yyaxis right
                     plot(trialData.mtime(2:end),vyTorso,'color',colors(2,:)),ylabel('Torso AP vel (m/s)');
                     set(gca,'ycolor',colors(2,:));
                     vline([trialData.mtime(indStart) trialData.mtime(indEnd)],'k-')
                     xlabel('Time (s)')
                     if plotind == 1 
                        titlename = sprintf('BL%i, Tr. %i %s',subj,trial,cond); 
                     else
                        titlename = sprintf('Tr. %i, %s',trial,cond); 
                     end
                     title(titlename);
                 end
            else
                 meanFz1(n) = nan;
                 meanFxy1(n) = nan;
                 meanFz2(n) = nan;
                 meanFxy2(n) = nan;
                 rho(n) = nan; pval(n) = nan;
            end
        end
        
        %% Concatenate metrics arrays depending on trial condition in this order of col's:  Decrease, SlowPost, PrefPre, Pref, PrefPost, Follow, Baseline, Increase, FastPost
        % Concatenate data only if trial is NOT bad
        if badTrial(n) ~= 1
            if strcmp(cond,'Decrease')
                speedDec = [speedDec speed(n)];
                LSLDec = [LSLDec LSLtrial(n)];
                RSLDec = [RSLDec RSL(n)];
                SLDec = [SLDec SL(n)];
                LcadDec = [LcadDec Lcad(n)];
                RcadDec = [RcadDec Rcad(n)];
                cadDec = [cadDec cad(n)];
                if analyzeForces == 1
                    meanFz1Dec = [meanFz1Dec meanFz1(n)];
                    meanFz2Dec = [meanFz2Dec meanFz2(n)];
                    corrFzDec = [corrFzDec rho(n)];
                    meanFxy1Dec = [meanFxy1Dec meanFxy1(n)];
                    meanFxy2Dec = [meanFxy2Dec meanFxy2(n)];
                end
            elseif strcmp(cond,'SlowPost')
                speedSlowPost = [speedSlowPost speed(n)];
                LSLSlowPost = [LSLSlowPost LSLtrial(n)];
                RSLSlowPost = [RSLSlowPost RSL(n)];
                SLSlowPost = [SLSlowPost SL(n)];
                LcadSlowPost = [LcadSlowPost Lcad(n)];
                RcadSlowPost = [RcadSlowPost Rcad(n)];
                cadSlowPost = [cadSlowPost cad(n)];
                if analyzeForces == 1
                    meanFz1SlowPost = [meanFz1SlowPost meanFz1(n)];
                    meanFz2SlowPost = [meanFz2SlowPost meanFz2(n)];
                    corrFzSlowPost = [corrFzSlowPost rho(n)];
                    meanFxy1SlowPost = [meanFxy1SlowPost meanFxy1(n)];
                    meanFxy2SlowPost = [meanFxy2SlowPost meanFxy2(n)];
                end
            elseif strcmp(cond,'Increase')
                speedInc = [speedInc speed(n)];
                LSLInc = [LSLInc LSLtrial(n)];
                RSLInc = [RSLInc RSL(n)];
                SLInc = [SLInc SL(n)];
                LcadInc = [LcadInc Lcad(n)];
                RcadInc = [RcadInc Rcad(n)];
                cadInc = [cadInc cad(n)];
                if analyzeForces == 1
                    meanFz1Inc = [meanFz1Inc meanFz1(n)];
                    meanFz2Inc = [meanFz2Inc meanFz2(n)];
                    corrFzInc = [corrFzInc rho(n)];
                    meanFxy1Inc = [meanFxy1Inc meanFxy1(n)];
                    meanFxy2Inc = [meanFxy2Inc meanFxy2(n)];
                end
            elseif strcmp(cond,'Preferred')
                speedPref = [speedPref speed(n)];
                LSLPref = [LSLPref LSLtrial(n)];
                RSLPref = [RSLPref RSL(n)];
                SLPref = [SLPref SL(n)];
                LcadPref = [LcadPref Lcad(n)];
                RcadPref = [RcadPref Rcad(n)];
                cadPref = [cadPref cad(n)];
                if analyzeForces == 1
                    meanFz1Pref = [meanFz1Pref meanFz1(n)];
                    meanFz2Pref = [meanFz2Pref meanFz2(n)];
                    corrFzPref = [corrFzPref rho(n)];
                    meanFxy1Pref = [meanFxy1Pref meanFxy1(n)];
                    meanFxy2Pref = [meanFxy2Pref meanFxy2(n)];
                end
            elseif strcmp(cond,'Follow')
                speedFollow = [speedFollow speed(n)];
                LSLFollow = [LSLFollow LSLtrial(n)];
                RSLFollow = [RSLFollow RSL(n)];
                SLFollow = [SLFollow SL(n)];
                LcadFollow = [LcadFollow Lcad(n)];
                RcadFollow = [RcadFollow Rcad(n)];
                cadFollow = [cadFollow cad(n)];
                if analyzeForces == 1
                    meanFz1Follow = [meanFz1Follow meanFz1(n)];
                    meanFz2Follow = [meanFz2Follow meanFz2(n)];
                    corrFzFollow = [corrFzFollow rho(n)];
                    meanFxy1Follow = [meanFxy1Follow meanFxy1(n)];
                    meanFxy2Follow = [meanFxy2Follow meanFxy2(n)];
                end
            elseif strcmp(cond,'PreferredPre')
                speedPrefPre = [speedPrefPre speed(n)];
                LSLPrefPre = [LSLPrefPre LSLtrial(n)];
                RSLPrefPre = [RSLPrefPre RSL(n)];
                SLPrefPre = [SLPrefPre SL(n)];
                LcadPrefPre = [LcadPrefPre Lcad(n)];
                RcadPrefPre = [RcadPrefPre Rcad(n)];
                cadPrefPre = [cadPrefPre cad(n)];
                if analyzeForces == 1
                    meanFz1PrefPre = [meanFz1PrefPre meanFz1(n)];
                    meanFz2PrefPre = [meanFz2PrefPre meanFz2(n)];
                    corrFzPrefPre = [corrFzPrefPre rho(n)];
                    meanFxy1PrefPre = [meanFxy1PrefPre meanFxy1(n)];
                    meanFxy2PrefPre = [meanFxy2PrefPre meanFxy2(n)];
                end
            elseif strcmp(cond,'PreferredPost')
                speedPrefPost = [speedPrefPost speed(n)];
                LSLPrefPost = [LSLPrefPost LSLtrial(n)];
                RSLPrefPost = [RSLPrefPost RSL(n)];
                SLPrefPost = [SLPrefPost SL(n)];
                LcadPrefPost = [LcadPrefPost Lcad(n)];
                RcadPrefPost = [RcadPrefPost Rcad(n)];
                cadPrefPost = [cadPrefPost cad(n)];
                if analyzeForces == 1
                    meanFz1PrefPost = [meanFz1PrefPost meanFz1(n)];
                    meanFz2PrefPost = [meanFz2PrefPost meanFz2(n)];
                    corrFzPrefPost = [corrFzPrefPost rho(n)];
                    meanFxy1PrefPost = [meanFxy1PrefPost meanFxy1(n)];
                    meanFxy2PrefPost = [meanFxy2PrefPost meanFxy2(n)];
                end
            elseif strcmp(cond,'FastPost')
                speedFastPost = [speedFastPost speed(n)];
                LSLFastPost = [LSLFastPost LSLtrial(n)];
                RSLFastPost = [RSLFastPost RSL(n)];
                SLFastPost = [SLFastPost SL(n)];
                LcadFastPost = [LcadFastPost Lcad(n)];
                RcadFastPost = [RcadFastPost Rcad(n)];
                cadFastPost = [cadFastPost cad(n)];
                if analyzeForces == 1
                    meanFz1FastPost = [meanFz1FastPost meanFz1(n)];
                    meanFz2FastPost = [meanFz2FastPost meanFz2(n)];
                    corrFzFastPost = [corrFzFastPost rho(n)];
                    meanFxy1FastPost = [meanFxy1FastPost meanFxy1(n)];
                    meanFxy2FastPost = [meanFxy2FastPost meanFxy2(n)];
                end

            elseif strcmp(cond,'Baseline')
                speedBase = [speedBase speed(n)];
                LSLBase = [LSLBase LSLtrial(n)];
                RSLBase = [RSLBase RSL(n)];
                SLBase = [SLBase SL(n)];
                LcadBase = [LcadBase Lcad(n)];
                RcadBase = [RcadBase Rcad(n)];
                cadBase = [cadBase cad(n)];
                if analyzeForces == 1
                    meanFz1Base = [meanFz1Base meanFz1(n)];
                    meanFz2Base = [meanFz2Base meanFz2(n)];
                    corrFzBase = [corrFzBase rho(n)];
                    meanFxy1Base = [meanFxy1Base meanFxy1(n)];
                    meanFxy2Base = [meanFxy2Base meanFxy2(n)];
                end
            end        
        end
    end
end

% %% Save metrics to mat file for subj
% if subj < 10
%     file = sprintf('BL0%i_processed',subj);
% else
%     file = sprintf('BL%i_processed',subj);
% end
% save file

%% Export and save fig(s)
if plotFz == 1
    % Format figure
    s2 = [5 1 8 6];
    condArray{1} = 'Inc';
    condArray{2} = 'Dec';
    condArray{3} = 'Pref';
    condArray{4} = 'Follow';
    for f = 1:4
        set(figure(f),'units', 'inches','paperunits','inches','pos',s2,'PaperOrientation','landscape');
        figname = sprintf('BL%i_Fz_%s',subj,condArray{f});
        saveas(figure(f),figname,'fig'); % Save as fig file
%         print(figname,'-dpdf','-fillpage')
%         export_fig(figname,'-pdf','-transparent','-append'); % Export to pdf
    end
end

%% Concatenate and save means and SDs for all trials of a cond combined in order of expected mean speed values:
% Decrease, SlowPost, PrefPre, Pref, PrefPost, Follow, Baseline (eyes open), Increase, FastPost
if subj < 4  
    mSpeed(1,:) = [nanmean(speedDec) nanstd(speedDec)];
    mSpeed(2,:) = [nanmean(speedSlowPost) nanstd(speedSlowPost)];
    mSpeed(3,:) = [nanmean(speedPrefPre) nanstd(speedPrefPre)];
    mSpeed(4,:) = [nanmean(speedPref) nanstd(speedPref)];
    mSpeed(5,:) = [nanmean(speedPrefPost) nanstd(speedPrefPost)];
    mSpeed(6,:) = [nanmean(speedFollow) nanstd(speedFollow)];
    mSpeed(7,:) = [nanmean(speedInc) nanstd(speedInc)];
    mSpeed(8,:) = [nanmean(speedFastPost) nanstd(speedFastPost)];
else
    mSpeed(1,:) = [nanmean(speedDec) nanstd(speedDec)];
    mSpeed(2,:) = [nanmean(speedSlowPost) nanstd(speedSlowPost)];
    mSpeed(3,:) = [nanmean(speedPrefPre) nanstd(speedPrefPre)];
    mSpeed(4,:) = [nanmean(speedPref) nanstd(speedPref)];
    mSpeed(5,:) = [nanmean(speedPrefPost) nanstd(speedPrefPost)];
    mSpeed(6,:) = [nanmean(speedFollow) nanstd(speedFollow)];
    mSpeed(7,:) = [nanmean(speedBase) nanstd(speedBase)];
    mSpeed(8,:) = [nanmean(speedInc) nanstd(speedInc)];
    mSpeed(9,:) = [nanmean(speedFastPost) nanstd(speedFastPost)];
end

mLSL(1,:) = [nanmean(LSLDec) nanstd(LSLDec)];
mLSL(2,:) = [nanmean(LSLPref) nanstd(LSLPref)];
mLSL(3,:) = [nanmean(LSLFollow) nanstd(LSLFollow)];
mLSL(4,:) = [nanmean(LSLInc) nanstd(LSLInc)];

mRSL(1,:) = [nanmean(RSLDec) nanstd(RSLDec)];
mRSL(2,:) = [nanmean(RSLPref) nanstd(RSLPref)];
mRSL(3,:) = [nanmean(RSLFollow) nanstd(RSLFollow)];
mRSL(4,:) = [nanmean(RSLInc) nanstd(RSLInc)];

mSL(1,:) = [nanmean(SLDec) nanstd(SLDec)];
mSL(2,:) = [nanmean(SLPref) nanstd(SLPref)];
mSL(3,:) = [nanmean(SLFollow) nanstd(SLFollow)];
mSL(4,:) = [nanmean(SLInc) nanstd(SLInc)];

mLcad(1,:) = [nanmean(LcadDec) nanstd(LcadDec)];
mLcad(2,:) = [nanmean(LcadPref) nanstd(LcadPref)];
mLcad(3,:) = [nanmean(LcadFollow) nanstd(LcadFollow)];
mLcad(4,:) = [nanmean(LcadInc) nanstd(LcadInc)];

mRcad(1,:) = [nanmean(RcadDec) nanstd(RcadDec)];
mRcad(2,:) = [nanmean(RcadPref) nanstd(RcadPref)];
mRcad(3,:) = [nanmean(RcadFollow) nanstd(RcadFollow)];
mRcad(4,:) = [nanmean(RcadInc) nanstd(RcadInc)];

mcad(1,:) = [nanmean(cadDec) nanstd(cadDec)];
mcad(2,:) = [nanmean(cadPref) nanstd(cadPref)];
mcad(3,:) = [nanmean(cadFollow) nanstd(cadFollow)];
mcad(4,:) = [nanmean(cadInc) nanstd(cadInc)];

% Also save cadences to test if patner pref matched solo prefj
mcadPref = [nanmean(cadPrefPre) nanmean(cadPref)];

filename = sprintf('BL%i_metrics.mat',subj);

% Save force metrics if they were valid (exclude BL04 and BL08)
if subj ~= 4 && subj ~= 6
    % Axial forces
    mFz1(1,:) = [nanmean(meanFz1Dec) nanstd(meanFz1Dec)];
    mFz1(2,:) = [nanmean(meanFz1Follow) nanstd(meanFz1Follow)];
    mFz1(3,:) = [nanmean(meanFz1Pref) nanstd(meanFz1Pref)];
    mFz1(4,:) = [nanmean(meanFz1Inc) nanstd(meanFz1Inc)];

    mFz2(1,:) = [nanmean(meanFz2Dec) nanstd(meanFz2Dec)];
    mFz2(2,:) = [nanmean(meanFz2Follow) nanstd(meanFz2Follow)];
    mFz2(3,:) = [nanmean(meanFz2Pref) nanstd(meanFz2Pref)];
    mFz2(4,:) = [nanmean(meanFz2Inc) nanstd(meanFz2Inc)];

    % Off-axis forces
    mFxy1(1,:) = [nanmean(meanFxy1Dec) nanstd(meanFxy1Dec)];
    mFxy1(2,:) = [nanmean(meanFxy1Follow) nanstd(meanFxy1Follow)];
    mFxy1(3,:) = [nanmean(meanFxy1Pref) nanstd(meanFxy1Pref)];
    mFxy1(4,:) = [nanmean(meanFxy1Inc) nanstd(meanFxy1Inc)];

    mFxy2(1,:) = [nanmean(meanFxy2Dec) nanstd(meanFxy2Dec)];
    mFxy2(2,:) = [nanmean(meanFxy2Follow) nanstd(meanFxy2Follow)];
    mFxy2(3,:) = [nanmean(meanFxy2Pref) nanstd(meanFxy2Pref)];
    mFxy2(4,:) = [nanmean(meanFxy2Inc) nanstd(meanFxy2Inc)];

    % Corr FT1 and FT2
    mCorrFz(1,:) = [nanmean(corrFzDec) nanstd(corrFzDec)];
    mCorrFz(2,:) = [nanmean(corrFzFollow) nanstd(corrFzFollow)];
    mCorrFz(3,:) = [nanmean(corrFzPref) nanstd(corrFzPref)];
    mCorrFz(4,:) = [nanmean(corrFzInc) nanstd(corrFzInc)];

    save(filename,'mSpeed','mSL','mcad','mFz1','mFz2','mFxy1','mFxy2','mCorrFz','biasFx1','biasFy1','biasFz1');
else
    save(filename,'mSpeed','mSL','mcad');
end
    
