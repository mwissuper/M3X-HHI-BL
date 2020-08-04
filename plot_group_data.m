% Load all participants and calc group metrics

clear; clc; close all;

subj_array = [2 3 5 7:11]; % only subj's with a zero F trial

for ind = 1:length(subj_array)
    subj = subj_array(ind);
    filename = sprintf('BL%i_metrics.mat',subj);
    data = load(filename);
    
    Fz1bias(ind) = data.biasFz1;
    xlab{ind} = num2str(subj);
end

plot(1:length(subj_array),Fz1bias,'x'),ylabel('Fz1 zero bias (V)');
box off, title('Solo walking trials as ref. for ea. participant');
set(gca,'tickdir','out','xticklabel',xlab),xlabel('Participant');
