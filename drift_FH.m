% Check drift on FH's when held horizontally

% Load trials
% Calculate mean voltage each force axis
% Plot volt vs. time and mean overlaid per trial to check constant voltage
% Plot mean volt. vs. trial to see drift for experiment. Maybe plot vs.
% time of trial too?

clear; clc; close all;

% Analysis options

subj = 11; 

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

%% Loop through all trials for participant
for n = 1:length(trialArray)
    trial = trialArray(n);
    filename = sprintf('calib0%i.mat',trialArray(n));
    trialData = load([sourcefolder,'\',subjFolder,'\',filename]);
     
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
    % change sampling to match marker data
    FT1.Fx = downsample(trialData.rawData.analog.other(:,chan.FT1.Fx),10); 
    FT1.Fy = downsample(trialData.rawData.analog.other(:,chan.FT1.Fy),10);
    FT1.Fz = downsample(trialData.rawData.analog.other(:,chan.FT1.Fz),10); 
    FT2.Fx = downsample(trialData.rawData.analog.other(:,chan.FT2.Fx),10); 
    FT2.Fy = downsample(trialData.rawData.analog.other(:,chan.FT2.Fy),10);
    FT2.Fz = downsample(trialData.rawData.analog.other(:,chan.FT2.Fz),10);
    
    indStart = 1; indEnd = 120; % Look at short time period and keep consistent across all trials
    
    % Median filter forces to reduce noise, multiply by factor to
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
    Fz1 = medfilt1(FT1.Fz)*SFz1;
    Fz2 = medfilt1(FT2.Fz)*SFz2;
    Fx1 = medfilt1(FT1.Fx)*SFxy1;
    Fx2 = medfilt1(FT2.Fx)*SFxy2;
    Fy1 = medfilt1(FT1.Fy)*SFxy1;
    Fy2 = medfilt1(FT2.Fy)*SFxy2;

    % Just use mean during period of interest as metric
    meanFz1(n) = nanmean(Fz1(indStart:indEnd));
    meanFx1(n) = nanmean(Fx1(indStart:indEnd));
    meanFy1(n) = nanmean(Fy1(indStart:indEnd));
    meanFz2(n) = nanmean(Fz2(indStart:indEnd));
    meanFx2(n) = nanmean(Fx2(indStart:indEnd));
    meanFy2(n) = nanmean(Fy2(indStart:indEnd));
    
    % Trial 6 has partner FT2 disconnected so do not use!
    if (subj == 10 && trial == 6) || (subj == 11 && trial == 5)
        meanFz2(n) = nan;
        meanFx2(n) = nan;
        meanFy2(n) = nan;
        Fz2 = nan; Fx2 = nan; Fy2 = nan;
    end

    %% Plots
    if plotFzTime == 1
        plotind = plotind + 1;
        subplot(numrows,numcols,plotind);
        if subj == 10 && trial ~= 6 || (subj == 11 && trial ~= 5)
            plot(trialData.mtime(indStart:indEnd),Fz1(indStart:indEnd),'g',trialData.mtime(indStart:indEnd),Fz2(indStart:indEnd),'b'); % Participant L hand blue, R hand green
        else
            plot(trialData.mtime(indStart:indEnd),Fz1(indStart:indEnd),'g'); % Participant FT1 only
        end
        hline(meanFz1(n),'g'); hline(meanFz2(n),'b');
        titlename = sprintf('Trial %i',trial);title(titlename);
        ylabel('Fz axial (N)');
    end
end

%% Plot mean vs. trial and linfit
[p1, stat] = polyfit(1:length(meanFz1),meanFz1,1);
subplot(2,1,1),plot(meanFz1,'gx'); ylabel('F (N)'),xlabel('Trial')
subplot(2,1,2),plot(meanFz2,'bx'); ylabel('F (N)'),xlabel('Trial')