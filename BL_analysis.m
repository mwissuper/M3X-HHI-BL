% Plot time series data for an individual subject. Calculate metrics for
% participant. Do all filtering in Nexus, so just look at
% rawData.video.markers for marker data

clear; clc; close all;

% Analysis options

subj = 3; % Total subj's 1:11. Do not trust FT1 for BL04 (looks not connected). BL01 force data might also be bad but need to check
analyzeForces = 0; % Tension is positive. Not confident of directions of off-axis force vectors unless use recoverForces so don't try to interpret directionality

%% Plotting options
plotTorso = 0; % 1 plot AP dir to check speed calculations
plotVyTorso = 1; % Check begninning of moving vs. not moving
plotHeel = 0;
plotFz = 0; % Plot force in N with HS events
plotFxy = 0;
plotFzBias = 0; % Plot Fz volt vs. time for each axis for bias voltage trial
plotFzV = 0; % Plot raw voltage for Fz for all trials as check on bias

if plotVyTorso == 1 || plotFzV == 1
    plotAllTrials = 1;
else
    plotAllTrials = 0; % plot only partner trials
end

pcount = 0;

%% Constants

if subj == 2
    trialArray = [1:49 51 53:56]; 
elseif subj == 3
    trialArray = [1:8 10:21 23:26 28 29 31:56 58 59]; % Some bad trials (22, 27, 30), trial 9 started with left foot and didn't catch it during expt. Incorrectly repeated pref in trial 57 instead of asymm
elseif subj == 5
    trialArray = 1:59;
else
    trialArray = 1:60;
end
sourcefolder = cd;
subjFolder = sprintf('BL0%i',subj);
 
configFile = sprintf('BL0%i_worksheet.mat',subj);
temp = load(configFile);
config = temp.data;

colors(1,:) = [0.00,0.45,0.74]; % nice blue
colors(2,:) = [0.85,0.33,0.10]; % nice red
colors(3,:) = [0.47,0.67,0.19]; % nice green

markerArray = {'CLAV','C7','LHEE','RHEE'}; % Markers to be analyzed

%% First get voltage bias for zero for FT1 for each participant. Use first solo walking trial
if analyzeForces == 1
    if subj == 2
        biasTrial = 45;
    elseif subj == 3
        biasTrial = 45; 
    elseif subj == 4 || subj == 5
        biasTrial = 49;
    end

    if biasTrial < 10
        filename = sprintf('Trial0%i.mat',biasTrial);
    else
        filename = sprintf('Trial%i.mat',biasTrial);
    end
    trialData = load([sourcefolder,'\',subjFolder,'\',filename]);
    s = size(trialData.rawData.analog.otherid);

    % Set up filter for force data
    fsA = trialData.AnalogFrameRate;
    wn = 60/(fsA/2);
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
    % change sampling to match marker data
    FT1.Fx = downsample(trialData.rawData.analog.other(:,chan.FT1.Fx),10); 
    FT1.Fy = downsample(trialData.rawData.analog.other(:,chan.FT1.Fy),10);
    FT1.Fz = downsample(trialData.rawData.analog.other(:,chan.FT1.Fz),10); 

    indStart = 1; indEnd = 120; % Look at short time period and keep consistent across all trials
    
    % Filter the force components
    Fx1 = filtfilt(Bf,Af,FT1.Fx1);
    Fy1 = filtfilt(Bf,Af,FT1.Fy1);
    Fz1 = filtfilt(Bf,Af,FT1.Fz1);

    % Just use mean during period of interest as metric
    biasFx1 = nanmean(Fx1(indStart:indEnd));
    biasFy1 = nanmean(Fy1(indStart:indEnd));
    biasFz1 = nanmean(Fz1(indStart:indEnd));
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
vyThresh = 0.001; % vel thresh to start trial (m/s)

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

for n = 14%1:length(trialArray)
    clear Markers
    trial = trialArray(n)
    if trial < 10
        filename = sprintf('Trial0%i.mat',trialArray(n));
    else
        filename = sprintf('Trial%i.mat',trialArray(n));
    end
    trialData = load([sourcefolder,'\',subjFolder,'\',filename]);
    
    fsM = trialData.VideoFrameRate;  
    
    % Set up filter for force data
    fsA = trialData.AnalogFrameRate;
    wn = 60/(fsA/2);
    [Bf,Af] = butter(4,wn,'low');
    
    % c3d files put in zeros for
    % gaps so need to find these and make these nan's!
    % Loop through all markers
    % Will only work for markers that are mostly labeled, so choose which
    % ones to preprocess this way.
    for i = 1:length(markerArray)
        indMarker(i) = findMarkerInd(markerArray{i},subj,trialData.MarkerID);
    end
        
    for i = indMarker
        Markers(:,i,:) = replaceZeros(trialData.rawData.video.markers(:,i,:));
%         if isempty(find(isnan(Markers(:,i,:)),1,'first')) % no nan values, ok to use filtfilt
%             Markers(:,i,:) = filtfilt(Bm,Am,Markers(:,i,:)); % will error out if there are nan's!
%         else
%             Markers(:,i,:) = filtIgnoreNan(Bm,Am,order,Markers(:,i,:));
%         end
    end

    cond = config{trialArray(n),2};
    if ~strcmp(cond,'Asymmetric')
        
        %% Find when forward walking started and ended based on clav marker
        
        % Find index of MarkerID for particular marker
        temp.indClav = findMarkerInd('CLAV',subj,trialData.MarkerID);
        clav = squeeze(Markers(:,temp.indClav,:));
        temp = [];
        clav = clav./1000; % convert to m
        
        % If gaps in clav, that means significant gap and can't use this
        % marker. Use C7 in this case.
        if ~isempty(find(isnan(clav),1,'first')) % sig. gaps
            temp.indC7 = findMarkerInd('C7',subj,trialData.MarkerID);
            c7 = squeeze(Markers(:,temp.indC7,:))./1000;
            temp = [];
            torso = c7;
            disp('c7')
        else
            torso = clav;
        end
        
        vyTorso = diff(torso(:,2)); 
        vyTorsoOffset = [vyTorso(2:end); nan];
        indStart = find(vyTorso > vyThresh,1,'first');
        indEnd = find(vyTorso(indStart:end) > vyThresh & vyTorsoOffset(indStart:end) < vyThresh,1,'last') + indStart;

        %% Plot vel of torso marker to check that beginning and end of fwd walking period found correctly
        if plotVyTorso == 1
            plotind = plotind + 1;
%             subplot(numcols,numrows,plotind)
            plot(trialData.mtime(2:end),vyTorso),ylabel('Torso AP vel (m/s)'); hold on;
            xlim([trialData.mtime(indStart)-0.2 trialData.mtime(indEnd)+0.2])
            vline([trialData.mtime(indStart) trialData.mtime(indEnd)],'k-')
            hline(vyThresh,'k--');
            titlename = sprintf('Trial %i',trial); title(titlename);
        end
        
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
            titlename = sprintf('Trial %i',trial); title(titlename);
        end
        
        %% Find HS events from heel marker data. Only use ones from forward walking period.
        % Find index of MarkerID for particular marker
       
        temp.indLheel = findMarkerInd('LHEE',subj,trialData.MarkerID);
        Lheel = squeeze(Markers(:,temp.indLheel,:));
%         temp = [];
        Lheel = Lheel./1000; % convert to m
        indLHS = getHS(Lheel(1:indEnd,3));
        indLHS(find(indLHS <= indStart | indLHS >= indEnd)) = [];
        
        temp.indRheel = findMarkerInd('RHEE',subj,trialData.MarkerID);
        Rheel = squeeze(Markers(:,temp.indRheel,:));
        temp = [];
        Rheel = Rheel./1000; % convert to m
        indRHS = getHS(Rheel(1:indEnd,3));
        indRHS(find(indRHS <= indStart | indRHS >= indEnd)) = [];
        
        % Don't look at the first initiation step(s) or last collection step
        % because not steady-state. Should be on R side.
        ind1 = find(indRHS >= indStart & indRHS <= indLHS(1));
        indRHS(ind1) = [];
        ind2 = find(indRHS >= indLHS(end) & indRHS <= indEnd);
        indRHS(ind2) = [];
        
        %% Commented out following bc done with old Nexus pipelines and old procBatch Marker output
%         if subj == 2 % exceptions where algo didn't work and manually correct event based on visual observation in Nexus
%             if trial == 18
%                 indRHS(2) = 561;
%             elseif trial == 32
%                 indLHS(3) = 1066;
%             end
%         elseif subj == 3
%             if trial == 5
%                 indRHS(1) = 408;
%                 indRHS(2) = 648;
%                 indLHS(1) = 281;
%             elseif trial == 15
%                 indRHS(1) = 381;
%             end
%         elseif subj == 4
%             if trial == 5
%                 indRHS(end) = 757;
%                 indLHS(end) = 846;
%             elseif trial == 9
%                 indRHS(1) = 514;
%             elseif trial == 14
%                 indLHS(end) = 710;
%             elseif trial == 16
%                 indLHS(end) = 764;
%             elseif trial == 24
%                 indRHS(2) = 721;
%                 indLHS(3) = 865;
%             elseif trial == 27
%                 indLHS(end-1) = [];
%             elseif trial == 31
%                 indLHS(3) = [];
%                 indRHS(1) = 425;
%                 indRHS(3) = 1082;
%             elseif trial == 34
%                 indLHS(3) = [];
%             elseif trial == 39
%                 indLHS(end) = 732;
%             elseif trial == 47
%                 indRHS(2) = 595;
%             end
%         elseif subj == 5
%             if trial == 40
%                 indRHS(2) = [];
%             end
%         end

        %% Calculate avg speed from distance traveled by torso marker during
        % forward walking portion
        dist(n) = sum(vyTorso(indStart:indEnd));
        speed(n) = dist(n)/(trialData.mtime(indEnd)-trialData.mtime(indStart));
        
        %% Plot marker pos to check HS found correctly
        if plotHeel == 1
            plotind = plotind + 1;
            subplot(numcols,numrows,plotind),
            hold on;
            plot(trialData.mtime,Lheel(:,3),'b-',trialData.mtime(indLHS),Lheel(indLHS,3),'bx');
            plot(trialData.mtime,Rheel(:,3),'g-',trialData.mtime(indRHS),Rheel(indRHS,3),'gx');
            ylabel('Heel vert pos (m)')
            vline([trialData.mtime(indStart) trialData.mtime(indEnd)],'k-')
            xlim([trialData.mtime(indStart)-0.2 trialData.mtime(indEnd)+0.2])
            titlename = sprintf('Trial %i',trial); title(titlename);
        end
        
        %% Calculate SL for L and R separately and combined
        
        LSLarray = calcSL(Lheel(:,2),Rheel(:,2),indLHS);
        RSLarray = calcSL(Rheel(:,2),Lheel(:,2),indRHS);
        SLarray = [LSLarray; RSLarray];
        LSLtrial(n) = nanmean(LSLarray);
        RSL(n) = nanmean(RSLarray);
        SL(n) = nanmean(SLarray);
        
        %% Calculate ST and cadence for L and R separately and combined
        
        LSTarray = calcST(indLHS,indRHS,fsM);
        RSTarray = calcST(indRHS,indLHS,fsM);
        STarray = [LSTarray RSTarray];
        LST(n) = nanmean(LSTarray);
        RST(n) = nanmean(RSTarray);
        ST(n) = nanmean(STarray);
        
        Lcad(n) = 1./LST(n).*60; 
        Rcad(n) = 1./RST(n).*60; 
        cad(n) = 1./ST(n).*60;
        
        %% Force calculations. 
        if analyzeForces == 1 
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
                FT1.Fx1 = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fx));
                FT1.Fy1 = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fy));
                FT1.Fz1 = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fz));
                
                % change sampling to match marker data
                Fx = downsample(FT1.Fx,fsa/fs); 
                Fy = downsample(FT1.Fy,fsa/fs);
                Fz = downsample(FT1.Fz,fsa/fs); 

                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                plot(trialData.mtime(indStart:indEnd),Fz1(indStart:indEnd),'g');
                hline(biasFz1,'g');
                if plotind == 1
                    titlename = sprintf('BL%i T%i %s',subj,n,config{n,2});
                else
                    titlename = sprintf('T%i %s',n,config{n,2});
                end
                ylim([-0.45 0.05]); box off; set(gca,'tickdir','out');
                title(titlename);
                ylabel('Fz1 (V)'); xlabel('Time (s)');
            end
            if ~strcmp(cond,'Baseline') & ~strcmp(cond,'PreferredPre') & ~strcmp(cond,'SlowPost') & ~strcmp(cond,'FastPost')
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
                 FT1.Fx1 = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fx));
                 FT1.Fy1 = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fy));
                 FT1.Fz1 = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fz));
                 FT2.Fx2 = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT2.Fx));
                 FT2.Fy2 = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT2.Fy));
                 FT2.Fz2 = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT2.Fz));
                
                 % change sampling to match marker data
                 Fx1 = downsample(FT1.Fx,fsa/fs); 
                 Fy1 = downsample(FT1.Fy,fsa/fs);
                 Fz1 = downsample(FT1.Fz,fsa/fs); 
                 Fx2 = downsample(FT2.Fx,fsa/fs); 
                 Fy2 = downsample(FT2.Fy,fsa/fs);
                 Fz2 = downsample(FT2.Fz,fsa/fs);

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
                 Fz1 = FT1.Fz*SFz1;
                 Fz2 = FT2.Fz*SFz2;
                 Fxy1 = Fxy1*SFxy1;
                 Fxy2 = Fxy2*SFxy2;
                 
                 % Just use mean during period of interest as metric
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
                        plot(trialData.mtime,Fz1,'color',colors(1,:)),hold on
                        plot(trialData.mtime,Fz2,'color',colors(2,:)); % Participant L hand FT2, R hand FT1
                        if plotind == 1
                            legend('R hand','L hand');
                        end
                        ylabel('F axial (N)');
                        if subj == 2
                            ylims = [-12 12];
                        elseif subj == 3
                            ylims = [-25 8];
                        elseif subj == 4
                            ylims = [-4 6];
                        elseif subj == 5
                            ylims = [-20 20];
                        end
                    elseif plotFxy == 1
                        plot(trialData.mtime,Fxy1,trialData.mtime,Fxy2); % Participant L hand FT2, R hand FT1
                        ylabel('F off-axis (N)');
                        if subj == 3
                            ylims = [-30 15];
                        elseif subj == 4
                            ylims = [-4 6];
                        end
                    end
%                     ylim(ylims);
                    hline(0,'k-');
                    if plotFz == 1 
                        h = hline(meanFz1(n),'--'); set(h,'color',colors(1,:));
                        h = hline(meanFz2(n),'--'); set(h,'color',colors(1,:)); 
                    elseif plotFxy == 1 
                        hline(meanFxy1(n),'b--');hline(meanFxy2(n),'r--'); 
                    end
                    xlim([trialData.mtime(indStart) trialData.mtime(indEnd)])
                    xlabel('Time (s)')
                    if plotind == 1 
                        titlename = sprintf('BL%i, Tr. %i %s, rho = %.2f',subj,trial,cond,rho(n)); 
                        legend('F1 R','F2 L');
                    else
                        titlename = sprintf('Tr. %i, rho = %.2f',trial,rho(n)); 
                    end
                    title(titlename);
                    v = vline(trialData.mtime(indLHS),'-'); set(v,'color',colors(1,:)); 
                    v = vline(trialData.mtime(indRHS),'-'); set(v,'color',colors(2,:)); % Participant L foot red, R foot blue
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

%% Export and save fig(s)
if plotFz == 1
    % Format figure
    s2 = [5 1 11 8.5];
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
%% This code is copied from old pipeline that processed all participants together
% but now it is being used for an individual subject. However, the code
% expects data to be in this format (all trials in one file), so
% concatenation is done here.
% if subj < 10
%     filename = sprintf('BL0%i_worksheet.mat',subj);
% else
%     filename = sprintf('BL%i_worksheet.mat',subj);
% end
% 
% if analyzeForces == 1 
%     % Get force data, format it so can use existing force/work code later
%     formatData = formatSingleTrial(trialData);
%     % Copy the formatted trial fields into the main experimental
%     % dataset
%     fields = fieldnames(formatData);
%     for k = 1:length(fields)
%         dataset(n).(fields{k}) = formatData.(fields{k});
%     end
%     allTrials = concatDataMat(filename,subjFolder); % Concatenate indiv trials to one file;
%     allTrialsPre = preprocessHHI2017(allTrials); % Reorganizes individual trial data
%     TrialData = mainWorkPowerAnalysisMW(allTrialsPre,subj); % Calculates work and power transfer
%     % Some don't have NOTES. Check for NOTES
%     for n = 1:length(TrialData)
%         if ~isfield(TrialData(n).Info,'Notes')
%             TrialData(n).Info.Notes = ' ';
%         end
%     end
%     filtTrialData = filterRepeatTrials(TrialData);
%     save(file,'TrialData');
% end

%% Concatenate and save means and SDs for all trials of a cond combined in order of expected mean speed values:
% Decrease, SlowPost, PrefPre, Pref, PrefPost, Follow, Baseline (eyes open), Increase, FastPost
if subj < 4  
    mSpeed(1,:) = [nanmean(speedDec) std(speedDec)];
    mSpeed(2,:) = [nanmean(speedSlowPost) std(speedSlowPost)];
    mSpeed(3,:) = [nanmean(speedPrefPre) std(speedPrefPre)];
    mSpeed(4,:) = [nanmean(speedPref) std(speedPref)];
    mSpeed(5,:) = [nanmean(speedPrefPost) std(speedPrefPost)];
    mSpeed(6,:) = [nanmean(speedFollow) std(speedFollow)];
    mSpeed(7,:) = [nanmean(speedInc) std(speedInc)];
    mSpeed(8,:) = [nanmean(speedFastPost) std(speedFastPost)];
else
    mSpeed(1,:) = [nanmean(speedDec) std(speedDec)];
    mSpeed(2,:) = [nanmean(speedSlowPost) std(speedSlowPost)];
    mSpeed(3,:) = [nanmean(speedPrefPre) std(speedPrefPre)];
    mSpeed(4,:) = [nanmean(speedPref) std(speedPref)];
    mSpeed(5,:) = [nanmean(speedPrefPost) std(speedPrefPost)];
    mSpeed(6,:) = [nanmean(speedFollow) std(speedFollow)];
    mSpeed(7,:) = [nanmean(speedBase) std(speedBase)];
    mSpeed(8,:) = [nanmean(speedInc) std(speedInc)];
    mSpeed(9,:) = [nanmean(speedFastPost) std(speedFastPost)];
end

mLSL(1,:) = [nanmean(LSLDec) std(LSLDec)];
mLSL(2,:) = [nanmean(LSLPref) std(LSLPref)];
mLSL(3,:) = [nanmean(LSLFollow) std(LSLFollow)];
mLSL(4,:) = [nanmean(LSLInc) std(LSLInc)];

mRSL(1,:) = [nanmean(RSLDec) std(RSLDec)];
mRSL(2,:) = [nanmean(RSLPref) std(RSLPref)];
mRSL(3,:) = [nanmean(RSLFollow) std(RSLFollow)];
mRSL(4,:) = [nanmean(RSLInc) std(RSLInc)];

mSL(1,:) = [nanmean(SLDec) std(SLDec)];
mSL(2,:) = [nanmean(SLPref) std(SLPref)];
mSL(3,:) = [nanmean(SLFollow) std(SLFollow)];
mSL(4,:) = [nanmean(SLInc) std(SLInc)];

mLcad(1,:) = [nanmean(LcadDec) std(LcadDec)];
mLcad(2,:) = [nanmean(LcadPref) std(LcadPref)];
mLcad(3,:) = [nanmean(LcadFollow) std(LcadFollow)];
mLcad(4,:) = [nanmean(LcadInc) std(LcadInc)];

mRcad(1,:) = [nanmean(RcadDec) std(RcadDec)];
mRcad(2,:) = [nanmean(RcadPref) std(RcadPref)];
mRcad(3,:) = [nanmean(RcadFollow) std(RcadFollow)];
mRcad(4,:) = [nanmean(RcadInc) std(RcadInc)];

mcad(1,:) = [nanmean(cadDec) std(cadDec)];
mcad(2,:) = [nanmean(cadPref) std(cadPref)];
mcad(3,:) = [nanmean(cadFollow) std(cadFollow)];
mcad(4,:) = [nanmean(cadInc) std(cadInc)];

% Also save cadences to test if patner pref matched solo prefj
mcadPref = [nanmean(cadPrefPre) nanmean(cadPref)];

% Axial forces
mFz1(1,:) = [nanmean(meanFz1Dec) std(meanFz1Dec)];
mFz1(2,:) = [nanmean(meanFz1Follow) std(meanFz1Follow)];
mFz1(3,:) = [nanmean(meanFz1Pref) std(meanFz1Pref)];
mFz1(4,:) = [nanmean(meanFz1Inc) std(meanFz1Inc)];

mFz2(1,:) = [nanmean(meanFz2Dec) std(meanFz2Dec)];
mFz2(2,:) = [nanmean(meanFz2Follow) std(meanFz2Follow)];
mFz2(3,:) = [nanmean(meanFz2Pref) std(meanFz2Pref)];
mFz2(4,:) = [nanmean(meanFz2Inc) std(meanFz2Inc)];

% Off-axis forces
mFxy1(1,:) = [nanmean(meanFxy1Dec) std(meanFxy1Dec)];
mFxy1(2,:) = [nanmean(meanFxy1Follow) std(meanFxy1Follow)];
mFxy1(3,:) = [nanmean(meanFxy1Pref) std(meanFxy1Pref)];
mFxy1(4,:) = [nanmean(meanFxy1Inc) std(meanFxy1Inc)];

mFxy2(1,:) = [nanmean(meanFxy2Dec) std(meanFxy2Dec)];
mFxy2(2,:) = [nanmean(meanFxy2Follow) std(meanFxy2Follow)];
mFxy2(3,:) = [nanmean(meanFxy2Pref) std(meanFxy2Pref)];
mFxy2(4,:) = [nanmean(meanFxy2Inc) std(meanFxy2Inc)];

% Corr FT1 and FT2
mCorrFz(1,:) = [nanmean(corrFzDec) std(corrFzDec)];
mCorrFz(2,:) = [nanmean(corrFzFollow) std(corrFzFollow)];
mCorrFz(3,:) = [nanmean(corrFzPref) std(corrFzPref)];
mCorrFz(4,:) = [nanmean(corrFzInc) std(corrFzInc)];

filename = sprintf('BL%i_processed.mat',subj);
save(filename,'mSpeed','mSL','mcad','mFz1','mFz2','mFxy1','mFxy2','mCorrFz');

%% Plot settings

if subj < 4 % didn't do baseline trials
    condName{1} = 'Pref Pre';
    condName{2} = 'Follow Partner';
    condName{3} = 'Pref Partner';
    condName{4} = 'Pref Post';
    condName{5} = 'Dec Partner';
    condName{6} = 'Slow Post';
    condName{7} = 'Inc Partner';
    condName{8} = 'Fast Post';
else
    condName{1} = 'Baseline';
    condName{2} = 'Pref Pre';
    condName{3} = 'Follow Partner';
    condName{4} = 'Pref Partner';
    condName{5} = 'Pref Post';
    condName{6} = 'Dec Partner';
    condName{7} = 'Slow Post';
    condName{8} = 'Inc Partner';
    condName{9} = 'Fast Post';
end

s = [992   898   761   513];

%% Plot speed metric
figure
errorbar(1:length(condName),mSpeed(:,1),mSpeed(:,2),'x')
set(gca,'xticklabel',condName); ylabel('Speed (m/s)');
xlim([0.5 length(condName)+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);
titlename = sprintf('BL%i',subj); title(titlename);

%% Plot settings
condName = [];
condName{1} = 'Decrease';
condName{2} = 'Follow';
condName{3} = 'Preferred';
condName{4} = 'Increase';
numConds = length(condName);

%% Plot SL metrics
% figure
figure
plotind = 0;

% plotind = plotind + 1;
% subplot(2,2,plotind)
% errorbar(1:numConds,mLSL(:,1),mLSL(:,2),'x')
% xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);
% set(gca,'xtick',1:numConds);
% set(gca,'xticklabel',condName); ylabel('LSL (m)');
% 
% plotind = plotind + 1;
% subplot(2,2,plotind)
% errorbar(1:numConds,mRSL(:,1),mRSL(:,2),'x')
% set(gca,'xticklabel',condName); ylabel('RSL (m)');
% set(gca,'xtick',1:numConds);
% xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);

plotind = plotind + 1;
subplot(2,2,plotind)
errorbar(1:numConds,mSL(:,1),mSL(:,2),'kx')
set(gca,'xticklabel',condName); ylabel('SL (m)');
set(gca,'xtick',1:numConds);
xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);
titlename = sprintf('BL%i',subj); title(titlename);

%% Plot cad metrics
% figure
% plotind = 0;

% plotind = plotind + 1;
% subplot(2,2,plotind)
% errorbar(1:numConds,mLcad(:,1),mLcad(:,2),'x')
% set(gca,'xticklabel',condName); ylabel('Lcad (steps/min)');
% set(gca,'xtick',1:numConds);
% xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);
% titlename = sprintf('BL%i',subj); title(titlename);
% 
% plotind = plotind + 1;
% subplot(2,2,plotind)
% errorbar(1:numConds,mRcad(:,1),mRcad(:,2),'x')
% set(gca,'xticklabel',condName); ylabel('Rcad (steps/min)');
% set(gca,'xtick',1:numConds);
% xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);

plotind = plotind + 1;
subplot(2,2,plotind)
errorbar(1:numConds,mcad(:,1),mcad(:,2),'kx')
set(gca,'xticklabel',condName); ylabel('cad (m)');
set(gca,'xtick',1:numConds);
xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);

%% Plot force metrics
% figure
% plotind = 0;

% plotind = plotind + 1;
% subplot(2,2,plotind)
% errorbar(1:numConds,mFz1(:,1),mFz1(:,2),'x')
% set(gca,'xticklabel',condName); ylabel('Fz1 R (N)');
% set(gca,'xtick',1:numConds);
% xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);
% 
% plotind = plotind + 1;
% subplot(2,2,plotind)
% errorbar(1:numConds,mFz2(:,1),mFz2(:,2),'x')
% set(gca,'xticklabel',condName); ylabel('Fz2 L (N)');
% set(gca,'xtick',1:numConds);
% xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);

% Overlay L and R Fz
plotind = plotind + 1;
subplot(2,2,plotind)
errorbar(1:numConds,mFz1(:,1),mFz1(:,2),'x'),hold on
errorbar(1:numConds,mFz2(:,1),mFz2(:,2),'x')
legend('R','L','orientation','horizontal');
set(gca,'xticklabel',condName); ylabel('Mean axial F (N)');
set(gca,'xtick',1:numConds);
xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);

% Overlay L and R Fxy
plotind = plotind + 1;
subplot(2,2,plotind)
errorbar(1:numConds,mFxy1(:,1),mFxy1(:,2),'x'),hold on
errorbar(1:numConds,mFxy2(:,1),mFxy2(:,2),'x')
set(gca,'xticklabel',condName); ylabel('Mean off-axis F (N)');
set(gca,'xtick',1:numConds);
xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);

