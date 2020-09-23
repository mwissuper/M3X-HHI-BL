% Load all participants and calc group metrics

clear; clc; close all;

subj_array = [2 3 5 7:11]; % only subj's with a zero F trial

for ind = 1:length(subj_array)
    subj = subj_array(ind);
    filename = sprintf('BL%i_metrics.mat',subj);
    data = load(filename);
    
    Fz1bias(ind) = data.biasFz1;
     % Convert to N after filtering. Multiply by factor to
     % convert to N from V. FT2 gains changed 12/5/19 from 5V
     % FS to 10V FS.
     if subj < 10
         SFz2 = 1/0.01;
     else
         SFz2 = 1/0.02;
     end
     SFz1 = 1/0.02;
     Fz1(ind) = data.biasFz1*SFz1; 
%      Fz2(ind) = data.biasFz2*SFz2;
                 
    xlab{ind} = num2str(subj);
end

% Plot in volts
plot(1:length(subj_array),Fz1bias,'x');
ylabel('Fz1 zero bias (V)');
box off, title('Solo walking trials as ref. for ea. participant');
set(gca,'tickdir','out','xticklabel',xlab),xlabel('Participant');

figure
% Plot in N
plot(1:length(subj_array),Fz1,'x');
ylabel('Fz1 zero bias (N)');
box off, title('Solo walking trials as ref. for ea. participant');
set(gca,'tickdir','out','xticklabel',xlab),xlabel('Participant');

% 1/4 of range of Fz drift across subj's is 1.43V. Use this for Fx or Fy in
% beamwalking study
