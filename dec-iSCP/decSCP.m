clc
clear all
close all

% Time settings and variables
T = 10; % Trajectory final time
h = 0.2; % time step duration
tk = 0:h:T;
K = T/h + 1; % number of time steps
Ts = 0.01; % period for interpolation @ 100Hz
t = 0:Ts:T; % interpolated time vector

% Initial positions
po1 = [-2,2,1.5];
po2 = [2,2,1.5];
po3 = [2,-2,1.5];
po4 = [-2,-2,1.5];
po5 = [-2,0,1.5];
po6 = [2,0,1.5];

po = cat(3,po1,po2,po3,po4,po5,po6);

% Final positions
pf1 = [2,-2,1.5];
pf2 = [-2,-2,1.5];
pf3 = [-2,2,1.5];
pf4 = [2,2,1.5];
pf5 = [2,0,1.5];
pf6 = [-2,0,1.5];

pf  = cat(3,pf1,pf2,pf3,pf4,pf5,pf6);

% Workspace boundaries
pmin = [-4,-4,0];
pmax = [4,4,2.5];

% Empty list of obstacles
l = [];

% Minimum distance between vehicles in m
rmin = 0.75;

% Maximum acceleration in m/s^2
alim = 1;

N = size(po,3); % number of vehicles

tic %measure the time it gets to solve the optimization problem
for i = 1:N 
    poi = po(:,:,i);
    pfi = pf(:,:,i);
    [pi, vi, ai] = singleiSCP(poi,pfi,h,K,pmin,pmax,rmin,alim,l);
    l = cat(3,l,pi);
    
    pk(:,:,i) = pi;
    vk(:,:,i) = vi;
    ak(:,:,i) = ai;
    
    % Interpolate solution with a 100Hz sampling
    p(:,:,i) = spline(tk,pi,t);
    v(:,:,i) = spline(tk,vi,t);
    a(:,:,i) = spline(tk,ai,t);
end
toc

%%
L = length(t);
colors = get(gca,'colororder');
colors = [colors; [1,0,0];[0,1,0];[0,0,1];[1,1,0];[0,1,1];...
           [0.5,0,0];[0,0.5,0];[0,0,0.5];[0.5,0.5,0]];
figure(1)
set(gcf,'currentchar',' ')
while get(gcf,'currentchar')==' '
    for k = 1:K
        for i = 1:N
            plot3(pk(1,k,i),pk(2,k,i),pk(3,k,i),'o', ...
                  'LineWidth',2, 'Color',colors(i,:));
            hold on;
            grid on;
            xlim([-4,4])
            ylim([-4,4])
            zlim([0,3.5])
            plot3(po(1,1,i), po(1,2,i), po(1,3,i),'^',...
                  'LineWidth',2,'Color',colors(i,:));
            plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i),'x',...
                  'LineWidth',2,'Color',colors(i,:));    
        end
        drawnow
    end
    pause(1)
    clf
end

%% Plotting
L = length(t);
for i = 1:N
    figure(1)
    plot3(p(1,:,i), p(2,:,i), p(3,:,i), 'LineWidth',1.5);
    hold on;
    grid on;
    xlim([-4,4])
    ylim([-4,4])
    zlim([0,3.5])
    plot3(po(1,1,i), po(1,2,i), po(1,3,i),'or','LineWidth',2);
    plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i),'xr','LineWidth',2);
    
    figure(2)
    diff = p(:,:,i) - repmat(pf(:,:,i),length(t),1)';
    dist = sqrt(sum(diff.^2,1));
    plot(t, dist, 'LineWidth',1.5);
    grid on;
    hold on;
    xlabel('t [s]')
    ylabel('Distance to target [m]');
    
    
    figure(3)
    subplot(3,1,1)
    plot(t,p(1,:,i),'LineWidth',1.5);
    plot(t,pmin(1)*ones(length(t),1),'--r','LineWidth',1.5);
    plot(t,pmax(1)*ones(length(t),1),'--r','LineWidth',1.5);
    ylabel('x [m]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,2)
    plot(t,p(2,:,i),'LineWidth',1.5);
    plot(t,pmin(2)*ones(length(t),1),'--r','LineWidth',1.5);
    plot(t,pmax(2)*ones(length(t),1),'--r','LineWidth',1.5);
    ylabel('y [m]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,3)
    plot(t,p(3,:,i),'LineWidth',1.5);
    plot(t,pmin(3)*ones(length(t),1),'--r','LineWidth',1.5);
    plot(t,pmax(3)*ones(length(t),1),'--r','LineWidth',1.5);
    ylabel('z [m]')
    xlabel ('t [s]')
    grid on;
    hold on;

    figure(4)
    subplot(3,1,1)
    plot(t,v(1,:,i),'LineWidth',1.5);
    ylabel('vx [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,2)
    plot(t,v(2,:,i),'LineWidth',1.5);
    ylabel('vy [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,3)
    plot(t,v(3,:,i),'LineWidth',1.5);
    ylabel('vz [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    figure(5)
    subplot(3,1,1)
    plot(t,a(1,:,i),'LineWidth',1.5);
    plot(t,alim*ones(length(t),1),'--r','LineWidth',1.5);
    plot(t,-alim*ones(length(t),1),'--r','LineWidth',1.5);
    ylabel('ax [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,2)
    plot(t,a(2,:,i),'LineWidth',1.5);
    plot(t,alim*ones(length(t),1),'--r','LineWidth',1.5);
    plot(t,-alim*ones(length(t),1),'--r','LineWidth',1.5);
    ylabel('ay [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,3)
    plot(t,a(3,:,i),'LineWidth',1.5);
    plot(t,alim*ones(length(t),1),'--r','LineWidth',1.5);
    plot(t,-alim*ones(length(t),1),'--r','LineWidth',1.5);
    ylabel('az [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;
   
end

figure(6)
for i = 1:N
    for j = 1:N
        if(i~=j)
            diff = p(:,:,i) - p(:,:,j);
            dist = sqrt(sum(diff.^2,1));
            plot(t, dist, 'LineWidth',1.5);
            grid on;
            hold on;
            xlabel('t [s]')
            ylabel('Inter-agent distance [m]');
        end
    end
end
plot(t,rmin*ones(length(t),1),'--r','LineWidth',1.5);