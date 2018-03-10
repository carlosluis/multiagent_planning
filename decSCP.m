clc
clear all
close all

T = 5;
h = 0.1;
tk = 0:h:5;
K = T/h + 1;
Ts = 0.01;
t = 0:Ts:T;

% Initial positions
po1 = [-1,0,0.5];
po2 = [0,0,0.5];
po = cat(3,po1,po2);

% Final positions
pf1 = [0,2,0.5];
pf2 = [-1,2,0.5];
pf  = cat(3, pf1, pf2);

% Workspace boundaries
pmin = [-4,-4,0];
pmax = [4,4,3.5];

% Empty 
l = [];

N = size(po,3);
tic
for i = 1:N 
%     diff =  - po(:,:,i);
    [pi, vi, ai] = singleSCP(po(:,:,i),pf(:,:,i),h,K,pmin,pmax,l);
    l = cat(3,l,pi);
    p(:,:,i) = spline(tk,pi,t);
    v(:,:,i) = spline(tk,vi,t);
    a(:,:,i) = spline(tk,ai,t);
end
toc

% Plotting
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
    subplot(3,1,1)
    plot(t,p(1,:,i),'LineWidth',1.5);
    ylabel('x [m]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,2)
    plot(t,p(2,:,i),'LineWidth',1.5);
    ylabel('y [m]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,3)
    plot(t,p(3,:,i),'LineWidth',1.5);
    ylabel('z [m]')
    xlabel ('t [s]')
    grid on;
    hold on;

    figure(3)
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

    figure(4)
    subplot(3,1,1)
    plot(t,a(1,:,i),'LineWidth',1.5);
    ylabel('ax [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,2)
    plot(t,a(2,:,i),'LineWidth',1.5);
    ylabel('ay [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,3)
    plot(t,a(3,:,i),'LineWidth',1.5);
    ylabel('az [m/s]')
    xlabel ('t [s]')
    grid on;
   
end

diff = p(:,:,1) - p(:,:,2);
dist = sqrt(sum(diff.^2,1));

figure(5)
plot(t, dist, 'LineWidth',1.5);
grid on;
xlabel('t [s]')
ylabel('distance [m]');




