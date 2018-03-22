clc
clear all
close all

% Time settings and variables
T = 12; % Trajectory final time
h = 0.2; % time step duration
tk = 0:h:T;
K = T/h + 1; % number of time steps
Ts = 0.01; % period for interpolation @ 100Hz
t = 0:Ts:T; % interpolated time vector
k_hor = 15;

% Initial positions
po1 = [-1,1,1.5];
po2 = [0,1,1.5];
po3 = [1,1,1.5];
po4 = [-1,0,1.5];
po5 = [0,0,1.5];
po6 = [1,0,1.5];
po7 = [-1,-1,1.5];
po8 = [0,-1,1.5];
po9 = [1,-1,1.5];

po = cat(3,po1,po2,po3,po4,po5,po6,po7,po8,po9);

% Final positions
pf1 = po8;
pf2 = po2;
pf3 = po4;
pf4 = po9;
pf5 = po5;
pf6 = po1;
pf7 = po6;
pf8 = po3;
pf9 = po7;

pf  = cat(3,pf1,pf2,pf3,pf4,pf5,pf6,pf7,pf8,pf9);

% Workspace boundaries
pmin = [-4,-4,0.2];
pmax = [4,4,2.5];

% Empty list of obstacles
l = [];

% Minimum distance between vehicles in m
rmin = 0.75;

% Maximum acceleration in m/s^2
alim = 1;

N = size(po,3); % number of vehicles

% Some pre computations
A = getPosMat(h,k_hor);
Aux = [1 0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];
A_initp = [];
A_init = eye(6);

Delta = getDeltaMat(k_hor); 

for k = 1:k_hor
    A_init = Aux*A_init;
    A_initp = [A_initp; A_init(1:3,:)];  
end

tic %measure the time it gets to solve the optimization problem
for k = 1:K
    for n = 1:N
        if k==1
            poi = po(:,:,n);
            pfi = pf(:,:,n);
            [pi,vi,ai] = initDMPC(poi,pfi,h,k_hor,K);
        else
            pok = pk(:,k-1,n);
            vok = vk(:,k-1,n);
            aok = ak(:,k-1,n);
            [pi,vi,ai] = solveDMPC(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta); 
        end
        new_l(:,:,n) = pi;
        pk(:,k,n) = pi(:,1);
        vk(:,k,n) = vi(:,1);
        ak(:,k,n) = ai(:,1);
    end
    l = new_l;
    pred(:,:,:,k) = l;
end      
       
toc

for i = 1:N
    p(:,:,i) = spline(tk,pk(:,:,i),t);
    v(:,:,i) = spline(tk,vk(:,:,i),t);
    a(:,:,i) = spline(tk,ak(:,:,i),t); 
end

%%
figure(1)
colors = get(gca,'colororder');
colors = [colors; [1,0,0];[0,1,0];[0,0,1];[1,1,0];[0,1,1];...
           [0.5,0,0];[0,0.5,0];[0,0,0.5];[0.5,0.5,0]];


set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'currentchar',' ')
while get(gcf,'currentchar')==' '
    for i = 1:N
    h_line(i) = animatedline('LineWidth',2,'Color',[0,0.8,0],'LineStyle','--');
    end
    for k = 1:K
        for i = 1:N
            clearpoints(h_line(i));
            addpoints(h_line(i),pred(1,:,i,k),pred(2,:,i,k),pred(3,:,i,k));
            drawnow
            hold on;
            grid on;
            xlim([-4,4])
            ylim([-4,4])
            zlim([0,3.5])
            plot3(pk(1,k,i),pk(2,k,i),pk(3,k,i),'o',...
                'LineWidth',2,'Color',colors(i,:));
            plot3(po(1,1,i), po(1,2,i), po(1,3,i),'^',...
                  'LineWidth',2,'Color',colors(i,:));
            plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i),'x',...
                  'LineWidth',2,'Color',colors(i,:));    
        end
%         pause(0.5)
    end
    clf
    pause(0.1)
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
