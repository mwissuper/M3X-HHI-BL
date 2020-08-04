% Check drift on FH's when held horizontally. See if bias values are
% comparable from solo walking trial vs. these trials. See if FT2 is
% consistent across BL10 and BL11.

% Load trials
% Calculate mean voltage each force axis
% Plot volt vs. time and mean overlaid per trial to check constant voltage
% Plot mean volt. vs. trial to see drift for experiment. Maybe plot vs.
% time of trial too?

clear; clc; close all;

%% Analysis options

subj = 11; 
fcA = 10;

%% Plotting options
plotFzTime = 1;
plotFxyTime = 0;

%% Constants

if subj == 10
    trialArray = 1:6;
elseif subj == 11
    trialArray = 1:5;
end
sourcefolder = cd;
subjFolder = sprintf('BL%i',subj);

numrows = 2;
numcols = 3;
plotind = 0;

indStart = 10; indEnd = 120; % Look at short time period and keep consistent across all trials. Shart transient at beg of trials

%% Loop through all trials for participant
for n = 1:length(trialArray)
    trial = trialArray(n);
    filename = sprintf('calib0%i.mat',trialArray(n));
    trialData = load([sourcefolder,'\',subjFolder,'\',filename]);
    
    % Set up filters
    fsM = trialData.VideoFrameRate;  
    fsA = trialData.AnalogFrameRate;
    wn = fcA/(fsA/2);
    [Bf,Af] = butter(4,wn,'low');
     
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
     FT1.Fx = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fx));
     FT1.Fy = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fy));
     FT1.Fz = filtfilt(Bf,Af,trialData.rawData.analog.other(:,chan.FT1.Fz));
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
    
    % Just use mean during period of interest as metric
    meanFz1(n) = nanmean(Fz1(indStart:indEnd));
    meanFz2(n) = nanmean(Fz2(indStart:indEnd));    
    
    % Trials where partner FT2 disconnected so do not use!
    if (subj == 10 && trial == 6) || (subj == 11 && trial == 5)
        meanFz2(n) = nan;
        meanFx2(n) = nan;
        meanFy2(n) = nan;
        Fz2 = nan; Fx2 = nan; Fy2 = nan;
    end

    %% Plots
    if plotFzTime == 1
        plotind = plotind + 1;
        subplot(numrows,numcols,plotind); hold on;
        if subj == 10 && trial ~= 6 || (subj == 11 && trial ~= 5)
            plot(trialData.mtime,Fz1);
            plot(trialData.mtime(indStart:indEnd),meanFz1(n)*ones(size(indStart:indEnd)),'linewidth',2);
            %,trialData.mtime(indStart:indEnd),Fz2(indStart:indEnd),'b'); % Participant L hand blue, R hand green
        else
            plot(trialData.mtime,Fz1);
            plot(trialData.mtime(indStart:indEnd),meanFz1(n)*ones(size(indStart:indEnd)),'linewidth',2); % Participant FT1 only
        end
        titlename = sprintf('Trial %i',trial);title(titlename);
        ylabel('Fz axial (V)');
        if subj == 10
            ylim([-0.07 -0.03]);
        end
    end
end

%% Plot mean vs. trial and mean across trials

titlename = sprintf('BL%i',subj);
if subj == 10
    xlab = {'8.5','18.5','28.5','38.5','51.5','63.5'};
    xtickarray = 1:6;
else
    xlab = {'18.5','28.5','38.5','49.5','60.5'};
    xtickarray = 1:5;
end

subplot(2,1,1),plot(meanFzV1,'g','marker','x'); ylabel('Fz1 (V)'),xlabel('Trial')
title(titlename); set(gca,'xtick',xtickarray,'xticklabel',xlab); box off; set(gca,'tickdir','out');
subplot(2,1,2),plot(meanFzV2,'b','marker','x'); ylabel('Fz2 (V)'),xlabel('Trial')
set(gca,'xtick',xtickarray,'xticklabel',xlab); box off; set(gca,'tickdir','out');
