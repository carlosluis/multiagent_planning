clc
close all
clear
load('comp_all_9.mat');
% green = [0.55,0.75,0.11];
green = [0.55,0.71,0]; %apple green
% Volumen of arena
V = (pmax(1)-pmin(1))*(pmax(2)-pmin(2))*(pmax(3)-pmin(3));
pi = 3.14159265359;
%%
% Probability of success plots
prob_dmpc = sum(success_dmpc,2)/trials*100;
prob_cup = sum(success_cup,2)/trials*100;
prob_dec = sum(success_dec,2)/trials*100;
figure(1)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([1:5]);
ylim([0,105])
xlim([0,5.5])
h1 = plot(N_vector/V,prob_cup,':b','Linewidth',2.5);
plot(N_vector/V,prob_cup,'ob', 'MarkerFaceColor', 'b','Linewidth',1,'markers',10);
h2 = plot(N_vector/V,prob_dec',':','Color',green,'Linewidth',2.5);
plot(N_vector/V,prob_dec,'o','Color',green, 'MarkerFaceColor', green,'Linewidth',1,'markers',10);
h3 = plot(N_vector/V,prob_dmpc,':r','Linewidth',2.5);
plot(N_vector/V,prob_dmpc,'or', 'MarkerFaceColor', 'r','Linewidth',1,'markers',10);
xlabel(' Workspace Density [agents/m³]')
ylabel('Success Probability [%]');
[h, icons, plots, s] = legend([h1,h2,h3],'Centralized','Decoupled','DMPC');
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');
%%
% Computation time
tmean_dmpc = nanmean(t_dmpc,2);
tstd_dmpc = nanstd(t_dmpc,1,2);
tmean_cup = nanmean(t_cup,2);
tstd_cup = nanstd(t_cup,1,2);
tmean_dec = nanmean(t_dec,2);
tstd_dec = nanstd(t_dec,1,2);
figure(2)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([1:5]);
ylim([0,1150])
xlim([0,5.5])
h1 = errorbar(N_vector/V, tmean_cup,tstd_cup,':b','LineWidth',2.5);
plot(N_vector/V,tmean_cup,'ob', 'MarkerFaceColor', 'b','Linewidth',1,'markers',10);
h2 = errorbar(N_vector/V, tmean_dec,tstd_dec,':','Color',green,'LineWidth',2.5);
plot(N_vector/V,tmean_dec,'o','Color',green, 'MarkerFaceColor', green,'Linewidth',1,'markers',10);
h3 = errorbar(N_vector/V, tmean_dmpc,tstd_dmpc,':r','LineWidth',2.5);
plot(N_vector/V,tmean_dmpc,'or', 'MarkerFaceColor', 'r','Linewidth',1,'markers',10);
message = sprintf('Successful \ntrials \nonly');
text(1,400,message,'Fontsize',20)
xlabel(' Workspace Density [agents/m³]')
ylabel('Computation time [s]');
[h, icons, plots, s] = legend([h1,h2,h3],'Centralized','Decoupled','DMPC');
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');
%%
% Average travelled distance
avg_dist_dmpc = nanmean(totdist_dmpc,2);
std_dist_dmpc = nanstd(totdist_dmpc,1,2);
avg_dist_cup = nanmean(totdist_cup,2);
std_dist_cup = nanstd(totdist_cup,1,2);
avg_dist_dec = nanmean(totdist_dec,2);
std_dist_dec = nanstd(totdist_dec,1,2);
figure(3)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([1:5]);
ylim([0,40])
xlim([0,5.5])
h1 = errorbar(N_vector/V, avg_dist_cup,std_dist_cup,':b','LineWidth', 2.5);

plot(N_vector/V,avg_dist_cup,'ob', 'MarkerFaceColor', 'b','Linewidth',1,'markers',10);
h2 = errorbar(N_vector/V, avg_dist_dec,std_dist_dec,':','Color',green,'LineWidth', 2.5);
plot(N_vector/V,avg_dist_dec,'o','Color',green,'MarkerFaceColor',green,'Linewidth',1,'markers',10);
h3 = errorbar(N_vector/V, avg_dist_dmpc,std_dist_dmpc,':r','LineWidth', 2.5);
plot(N_vector/V,avg_dist_dmpc,'or', 'MarkerFaceColor', 'r','Linewidth',1,'markers',10);
text(0.5,18,message,'Fontsize',20)
xlabel(' Workspace Density [agents/m³]')
ylabel('Total Distance [m]');
[h, icons, plots, s] = legend([h1,h2,h3],'Centralized','Decoupled','DMPC');
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
h1 = errorbar(N_vector/V,tmean_traj,tstd_traj,':r','LineWidth', 2.5);
plot(N_vector/V,tmean_traj,'or', 'MarkerFaceColor', 'r','Linewidth',1.5,'markers',10);
text(0.5,10,message,'Fontsize',20)
xlabel(' Workspace Density [agents/m³]')
ylabel('Transition Time [s]');
[h, icons, plots, s] = legend([h1],'DMPC');
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');
