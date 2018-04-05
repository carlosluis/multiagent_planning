clc
clear
close all
load('26_comp.mat')
% Post-Processing
% Probability of success plots
prob_dec = sum(success_dec,2)/trials;
prob_dmpc = sum(success_dmpc,2)/trials;
figure(1)
plot(N_vector,prob_dec','r','Linewidth',2);
box on;
hold on;
ylim([0,1.05])
xlim([0,27])
plot(N_vector,prob_dmpc,'b','Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Success Probability');
legend('dec-iSCP','DMPC');

% Computation time
tmean_dec = nanmean(t_dec,2);
tstd_dec = nanstd(t_dec,1,2);
tmean_dmpc = nanmean(t_dmpc,2);
tstd_dmpc = nanstd(t_dmpc,1,2);
figure(2)
errorbar(N_vector,tmean_dec,tstd_dec,'r','Linewidth',2);
box on;
hold on;
xlim([0,27])
errorbar(N_vector,tmean_dmpc,tstd_dmpc,'b','Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Average Computation Time [s]');
legend('dec-iSCP','DMPC');