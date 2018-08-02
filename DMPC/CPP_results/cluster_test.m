clc
close all
clear
M = dlmread('cluster_test.txt','');
clus_arr_size = M(1,1);
num_veh_size = M(1,2);
num_trials = M(1,3);

cluster_size = M(2,1:clus_arr_size);
num_vehicles = M(2,clus_arr_size + 1 : clus_arr_size + num_veh_size);

for i=1:clus_arr_size
    for j = 1:num_veh_size
        time(i,j,:) = M(3 + i*j - 1, 1:num_trials);
    end
end

for i = 1:clus_arr_size
    avg_time(i) = mean(squeeze(time(i,:,:)),1);
    std_time(i) = std(squeeze(time(i,:,:)),1,1);
end

%% Plotting
colors = distinguishable_colors(clus_arr_size);

for i = 1:clus_arr_size
    figure(1)
    h_plot(i) = errorbar(num_vehicles,avg_time(:,i),std_time(:,i), 'LineWidth',1.5,...
                'Color',colors(i,:));
    h_label{i} = [num2str(cluster_size(i)) 'cluster(s)'];
    hold on;
    grid on;
    plot(num_vehicles, avg_time(:,i),'o', 'MarkerFaceColor', colors(i,:),...
        'Linewidth',1,'Color',colors(i,:));
    xlabel('Number of Vehicles');
    ylabel('Average Computation Time [s]');  
end
legend(h_plot,h_label);