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
k_hor = 10;

% Initial positions
po1 = [-2,0,1.5];
po2 = [2,0,1.5];

po = cat(3,po1,po2);

% Final positions
pf1 = [2,0,1.5];
pf2 = [-2,0,1.5];

pf  = cat(3,pf1,pf2);

% Workspace boundaries
pmin = [-4,-4,0];
pmax = [4,4,2.5];

% Empty list of obstacles
l = [];

% Minimum distance between vehicles in m
rmin = 0.9;

% Maximum acceleration in m/s^2
alim = 1;

N = size(po,3); % number of vehicles

tic %measure the time it gets to solve the optimization problem

for k = 1:K
    for n = 1:N
        if k==1
            poi = po(:,:,n);
            pfi = pf(:,:,n);
            pi = initSolution(poi,pfi,h,K);
            pi = pi(:,1:k_hor);
            vi = zeros(3,K);
            ai = zeros(3,K);
            l(:,:,n) = pi;
            pk(:,k,n) = pi(:,1);
            vk(:,k,n) = vi(:,1);
            ak(:,k,n) = ai(:,1); 
        else
            pok = pk(:,k-1,n);
            vok = vk(:,k-1,n);
            aok = ak(:,k-1,n);
            [pi,vi,ai] = solveDMPC(pok',pf(:,:,n),vok',aok',n,h,k,l,k_hor,rmin,pmin,pmax,alim); 
            l(:,:,n) = pi;
            pk(:,k,n) = pi(:,2);
            vk(:,k,n) = vi(:,2);
            ak(:,k,n) = ai(:,2); 
        end
           
    end
end
            
       
toc

L = length(t);
colors = get(gca,'colororder');
colors = [colors; [1,0,0];[0,1,0];[0,0,1]];
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

figure(1)
for i = 1:N 
    plot3(p(1,:,i), p(2,:,i), p(3,:,i), 'LineWidth',1.5);    
    hold on;
    grid on;
    xlim([-4,4])
    ylim([-4,4])
    zlim([0,3.5])
    plot3(po(1,1,i), po(1,2,i), po(1,3,i),'or','LineWidth',2);
    plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i),'xr','LineWidth',2);
end


% Plotting
% L = length(t);
% for i = 1:N
%     figure(1)
%     plot3(p(1,:,i), p(2,:,i), p(3,:,i), 'LineWidth',1.5);
%     hold on;
%     grid on;
%     xlim([-4,4])
%     ylim([-4,4])
%     zlim([0,3.5])
%     plot3(po(1,1,i), po(1,2,i), po(1,3,i),'or','LineWidth',2);
%     plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i),'xr','LineWidth',2); 
%     
%     figure(2)
%     subplot(3,1,1)
%     plot(t,p(1,:,i),'LineWidth',1.5);
%     ylabel('x [m]')
%     xlabel ('t [s]')
%     grid on;
%     hold on;
% 
%     subplot(3,1,2)
%     plot(t,p(2,:,i),'LineWidth',1.5);
%     ylabel('y [m]')
%     xlabel ('t [s]')
%     grid on;
%     hold on;
% 
%     subplot(3,1,3)
%     plot(t,p(3,:,i),'LineWidth',1.5);
%     ylabel('z [m]')
%     xlabel ('t [s]')
%     grid on;
%     hold on;
% 
%     figure(3)
%     subplot(3,1,1)
%     plot(t,v(1,:,i),'LineWidth',1.5);
%     ylabel('vx [m/s]')
%     xlabel ('t [s]')
%     grid on;
%     hold on;
% 
%     subplot(3,1,2)
%     plot(t,v(2,:,i),'LineWidth',1.5);
%     ylabel('vy [m/s]')
%     xlabel ('t [s]')
%     grid on;
%     hold on;
% 
%     subplot(3,1,3)
%     plot(t,v(3,:,i),'LineWidth',1.5);
%     ylabel('vz [m/s]')
%     xlabel ('t [s]')
%     grid on;
%     hold on;
% 
%     figure(4)
%     subplot(3,1,1)
%     plot(t,a(1,:,i),'LineWidth',1.5);
%     ylabel('ax [m/s]')
%     xlabel ('t [s]')
%     grid on;
%     hold on;
% 
%     subplot(3,1,2)
%     plot(t,a(2,:,i),'LineWidth',1.5);
%     ylabel('ay [m/s]')
%     xlabel ('t [s]')
%     grid on;
%     hold on;
% 
%     subplot(3,1,3)
%     plot(t,a(3,:,i),'LineWidth',1.5);
%     ylabel('az [m/s]')
%     xlabel ('t [s]')
%     grid on;
%    
% end
% 
% diff = p(:,:,6) - p(:,:,10);
% dist = sqrt(sum(diff.^2,1));
% 
% figure(5)
% plot(t, dist, 'LineWidth',1.5);
% grid on;
% xlabel('t [s]')
% ylabel('distance [m]');