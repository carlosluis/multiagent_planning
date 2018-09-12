clc
close all
clear
M = dlmread('cluster_test(100-ooqp).txt','');
clus_arr_size = M(1,1);
num_veh_size = M(1,2);
num_trials = M(1,3);

cluster_size = M(2,1:clus_arr_size);
num_vehicles = M(2,clus_arr_size + 1 : clus_arr_size + num_veh_size);

for i=1:clus_arr_size
    for j = 1:num_veh_size
        time(i,j,:) = M(3 + (i-1)*num_veh_size + j-1, 1:num_trials);
    end
end

for i = 1:clus_arr_size
    avg_time(i,:) = mean(squeeze(time(i,:,:)),2);
    std_time(i,:) = std(squeeze(time(i,:,:)),1,2);
end
clusters2plot = [1,2,4,8];
%% Plotting
colors = distinguishable_colors(clus_arr_size);
colors(3,:) = [0.55,0.75,0.11];
idx = 1;
for i = 1:clus_arr_size
    figure(1)
    xlim([0 110])
    set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
    set(gca,'FontSize',18)
    if (ismember(cluster_size(i),clusters2plot))
        h_plot(idx) = errorbar(num_vehicles,avg_time(i,:),std_time(i,:),'--', 'LineWidth',1.5,...
                    'Color',colors(idx,:));
        h_label{idx} = [num2str(cluster_size(i)) ' cluster(s)'];
        hold on;
        box on;
        plot(num_vehicles, avg_time(i,:),'o', 'MarkerFaceColor', colors(idx,:),...
            'Linewidth',1,'Color',colors(idx,:));
        xlabel('Number of Agents');
        ylabel('Computation Time [s]');
        idx = idx + 1;
    end
end
[h, icons, plots, s] = legend(h_plot,h_label);
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');