clc
close all
clear
load('comp_hard_soft_OnDemand_200.mat');
% green = [0.55,0.75,0.11];
green = [0.55,0.71,0]; %apple green
% Volumen of arena
V = (pmax(1)-pmin(1))*(pmax(2)-pmin(2))*(pmax(3)-pmin(3));
pi = 3.14159265359;
%%
% Probability of success plots
prob_dmpc = sum(success_dmpc,2)/trials*100;
prob_dmpc2 = sum(success_dmpc2,2)/trials*100;
prob_dmpc3 = sum(success_dmpc3,2)/trials*100;
figure(1)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([20 60 100 150 200]);
ylim([0,110])
xlim([0,220])
h1 = plot(N_vector,prob_dmpc2,':b','Linewidth',2.5);
plot(N_vector,prob_dmpc2,'ob', 'MarkerFaceColor', 'b','Linewidth',1,'markers',10);
h2 = plot(N_vector,prob_dmpc3,':','Color',green,'Linewidth',2.5);
plot(N_vector,prob_dmpc3,'o','Color',green, 'MarkerFaceColor', green,'Linewidth',1,'markers',10);
h3 = plot(N_vector,prob_dmpc,':r','Linewidth',2.5);
plot(N_vector,prob_dmpc,'or', 'MarkerFaceColor', 'r','Linewidth',1,'markers',10);
xlabel('Number of agents')
ylabel('Success Probability [%]');
[h, icons, plots, s] = legend([h1,h2,h3],'Hard','Hard On-Demand','Soft On-Demand');
set(h,'color','none');
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');
%%
% Computation time
tmean_dmpc = nanmean(t_dmpc,2);
tstd_dmpc = nanstd(t_dmpc,1,2);
tmean_dmpc2 = nanmean(t_dmpc2,2);
tstd_dmpc2 = nanstd(t_dmpc2,1,2);
tmean_dmpc3 = nanmean(t_dmpc3,2);
tstd_dmpc3 = nanstd(t_dmpc3,1,2);
figure(2)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([20 60 100 150 200]);
ylim([0,320])
xlim([0,220])
h1 = errorbar(N_vector, tmean_dmpc2,tstd_dmpc2,':b','LineWidth',2.5);
plot(N_vector,tmean_dmpc2,'ob', 'MarkerFaceColor', 'b','Linewidth',1,'markers',10);
h2 = errorbar(N_vector, tmean_dmpc3,tstd_dmpc3,':','Color',green,'LineWidth',2.5);
plot(N_vector,tmean_dmpc3,'o','Color',green, 'MarkerFaceColor', green,'Linewidth',1,'markers',10);
h3 = errorbar(N_vector, tmean_dmpc,tstd_dmpc,':r','LineWidth',2.5);
plot(N_vector,tmean_dmpc,'or', 'MarkerFaceColor', 'r','Linewidth',1,'markers',10);
message = sprintf('Successful \ntrials \nonly');
text(10,100,message,'Fontsize',20)
xlabel('Number of agents')
ylabel('Computation time [s]');
[h, icons, plots, s] = legend([h1,h2,h3],'Hard','Hard On-Demand','Soft On-Demand');
set(h,'color','none');
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');

%%
% Trajectory time
tmean_traj = nanmean(traj_time,2);
tstd_traj = nanstd(traj_time,1,2);
figure(4)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([1:5]);
xlim([0,6])
ylim([0,15])
h1 = errorbar(N_vector/V,tmean_traj,tstd_traj,':r','LineWidth', 2.5);
plot(N_vector/V,tmean_traj,'or', 'MarkerFaceColor', 'r','Linewidth',1.5,'markers',10);
text(0.5,8,message,'Fontsize',20)
xlabel(' Workspace Density [agents/m³]')
ylabel('Transition Time [s]');
[h, icons, plots, s] = legend([h1],'DMPC');
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');
