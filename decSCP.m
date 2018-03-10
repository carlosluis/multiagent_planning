clc
clear all
close all

T = 5;
h = 0.2;
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
pf2 = [1,2,0.5];
pf  = cat(3, pf1, pf2);

% Workspace boundaries
pmin = [-4,-4,0];
pmax = [4,4,3.5];

% Empty 
l = [];

N = size(po,3);

for i = 1:N 
%     diff =  - po(:,:,i);
    [pi, vi, ai] = singleSCP(po(:,:,i),pf(:,:,i),h,K,pmin,pmax,l);
    l = cat(3,l,pi);
    p(:,:,i) = spline(tk,pi,t);
    v(:,:,i) = spline(tk,vi,t);
    a(:,:,i) = spline(tk,ai,t);
end

% Plotting
L = length(t);
for i = 1:N
    figure(1)
    plot3(p(1,:,i), p(2,:,i), p(3,:,i));
    hold on;
    grid on;
    plot3(p(1,1,i), p(2,1,i), p(3,1,i),'or','LineWidth',2);
    plot3(p(1,L,i), p(2,L,i), p(3,L,i),'xr','LineWidth',2); 
    
    figure(2)
    subplot(3,1,1)
    plot(t,p(1,:,i));
    ylabel('x [m]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,2)
    plot(t,p(2,:,i));
    ylabel('y [m]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,3)
    plot(t,p(3,:,i));
    ylabel('z [m]')
    xlabel ('t [s]')
    grid on;
    hold on;

    figure(3)
    subplot(3,1,1)
    plot(t,v(1,:,i));
    ylabel('vx [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,2)
    plot(t,v(2,:,i));
    ylabel('vy [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,3)
    plot(t,v(3,:,i));
    ylabel('vz [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    figure(4)
    subplot(3,1,1)
    plot(t,a(1,:,i));
    ylabel('ax [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,2)
    plot(t,a(2,:,i));
    ylabel('ay [m/s]')
    xlabel ('t [s]')
    grid on;
    hold on;

    subplot(3,1,3)
    plot(t,a(3,:,i));
    ylabel('az [m/s]')
    xlabel ('t [s]')
    grid on;
   
end




