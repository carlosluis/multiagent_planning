clc
close all
clear
load('failure_rate3.mat');

% Failure analysis
violation_num = sum(violation,2);
goal_num = sum(failed_goal,2);
infes_num = sum(~feasible,2);
total_num = sum(violation_num) + sum(goal_num) + sum(infes_num);

StackData = [infes_num/trials*100 violation_num/trials*100 goal_num/trials*100];

figure(1)
h=bar(N_vector,StackData,'stacked','BarWidth',0.3);
myC= summer(size(StackData,2));
% for k = 1:size(StackData,2)
%     set(h(k),'facecolor',myC(k,:));
% end
box on;

xlim([0 210])
set(gca,'FontSize',18)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'XTickLabel',{[20] "" "" [80] "" "" [140] "" "" [200]})
xlabel('Number of agents');
ylabel('Failure Rate [%]');
legend('Collision during execution','Collision after interpolation')