clc
close all
clear
M = dlmread('cluster_test.txt','');
clus_arr_size = M(1,1);
num_veh_size = M(1,2);

cluster_size = M(2,1:clus_arr_size);
num_vehicles = M(2,clus_arr_size + 1 : clus_arr_size + num_veh_size);

avg_time = M(3:2+clus_arr_size, 1: num_veh_size);
%% Plotting
colors = distinguishable_colors(clus_arr_size);

for i = 1:clus_arr_size
    figure(1)
    h_plot(i) = plot(num_vehicles, avg_time(i,:), 'LineWidth',1.5,...
                'Color',colors(i,:));
    h_label{i} = [num2str(cluster_size(i)) 'cluster(s)'];
    hold on;
    grid on;
    xlabel('Number of Vehicles');
    ylabel('Average Computation Time [s]');  
end
legend(h_plot,h_label);