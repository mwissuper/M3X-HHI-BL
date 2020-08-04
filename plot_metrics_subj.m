% After running BL_analysis, plot summary metrics per subj

clear; clc; close all;

subj = 11; 

filename = sprintf('BL%i_metrics.mat',subj);
load(filename);

%% Plot settings

% Decrease, SlowPost, PrefPre, Pref, PrefPost, Follow, Baseline (eyes open), Increase, FastPost
if subj < 4 % didn't do baseline trials
    condName{1} = 'Dec Partner';
    condName{2} = 'Slow Post';
    condName{3} = 'Pref Pre';
    condName{4} = 'Pref Partner';
    condName{5} = 'Pref Post';
    condName{6} = 'Follow Partner';
    condName{7} = 'Inc Partner';
    condName{8} = 'Fast Post';
    indConds = [1 6 4 7];% For plots of partner cond's
else
    condName{1} = 'Dec Partner';
    condName{2} = 'Slow Post';
    condName{3} = 'Pref Pre';
    condName{4} = 'Pref Partner';
    condName{5} = 'Pref Post';
    condName{6} = 'Follow Partner';
    condName{7} = 'Baseline';
    condName{8} = 'Inc Partner';
    condName{9} = 'Fast Post';
    indConds = [1 6 4 8];% For plots of partner cond's
end

s = [992   898   761   513];

%% Plot speed metric (all speeds for all cond's)
% figure
% errorbar(1:length(condName),mSpeed(:,1),mSpeed(:,2),'x')
% set(gca,'xticklabel',condName); ylabel('Speed (m/s)');
% xlim([0.5 length(condName)+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);
% titlename = sprintf('BL%i',subj); title(titlename);

%% Plot settings
condName = [];
condName{1} = 'Decrease';
condName{2} = 'Follow';
condName{3} = 'Preferred';
condName{4} = 'Increase';
numConds = length(condName);

%% Plot speed metric
figure
plotind = 0;

plotind = plotind + 1;
subplot(2,2,plotind)
errorbar(1:numConds,mSpeed(indConds,1),mSpeed(indConds,2),'kx')
set(gca,'xticklabel',condName); ylabel('Speed (m/s)');
xlim([0.5 length(condName)+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);
titlename = sprintf('BL%i',subj); title(titlename);

%% Plot SL metrics
% figure

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
close all; figure
% plotind = 0;
if subj ~= 4 && subj ~= 6
    plotind = plotind + 1;
%     subplot(2,2,plotind)
    errorbar(1:numConds,mFz1(:,1),mFz1(:,2),'kx')
    set(gca,'xticklabel',condName); ylabel('Fz1 R (N)');
    set(gca,'xtick',1:numConds);
    xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);
    hline(0,'k--');
end
titlename = sprintf('BL%i',subj); title(titlename);
set(gcf,'outerposition',[541   514   300   284]);

% s = [108   499   597   457];
% set(gcf,'outerposition',s);
% 
% plotind = plotind + 1;
% subplot(2,2,plotind)
% errorbar(1:numConds,mFz2(:,1),mFz2(:,2),'x')
% set(gca,'xticklabel',condName); ylabel('Fz2 L (N)');
% set(gca,'xtick',1:numConds);
% xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);

% % Overlay L and R Fz
% plotind = plotind + 1;
% subplot(2,2,plotind)
% errorbar(1:numConds,mFz1(:,1),mFz1(:,2),'x'),hold on
% errorbar(1:numConds,mFz2(:,1),mFz2(:,2),'x')
% legend('R','L','orientation','horizontal');
% set(gca,'xticklabel',condName); ylabel('Mean axial F (N)');
% set(gca,'xtick',1:numConds);
% xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);
% 
% % Overlay L and R Fxy
% plotind = plotind + 1;
% subplot(2,2,plotind)
% errorbar(1:numConds,mFxy1(:,1),mFxy1(:,2),'x'),hold on
% errorbar(1:numConds,mFxy2(:,1),mFxy2(:,2),'x')
% set(gca,'xticklabel',condName); ylabel('Mean off-axis F (N)');
% set(gca,'xtick',1:numConds);
% xlim([0.5 numConds+0.5]); box off; set(gca,'tickdir','out'), set(gcf,'outerposition',s);

