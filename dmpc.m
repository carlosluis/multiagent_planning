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
k_hor = 5;

% Initial positions
po1 = [-2,0.01,1.5];
po2 = [2,0,1.5];

po = cat(3,po1,po2);

% Final positions
pf1 = [2,0.01,1.5];
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
            new_l(:,:,n) = pi;
            pk(:,k,n) = pi(:,1);
            vk(:,k,n) = vi(:,1);
            ak(:,k,n) = ai(:,1); 
        else
            pok = pk(:,k-1,n);
            vok = vk(:,k-1,n);
            aok = ak(:,k-1,n);
            [pi,vi,ai] = solveDMPC(pok',pf(:,:,n),vok',aok',n,h,k,l,k_hor,rmin,pmin,pmax,alim); 
            new_l(:,:,n) = pi;
            pk(:,k,n) = pi(:,2);
            vk(:,k,n) = vi(:,2);
            ak(:,k,n) = ai(:,2);
        end
    end
    l = new_l;
    pred(:,:,:,k) = l;
end      
       
toc
%%
for i = 1:N
    p(:,:,i) = spline(tk,pk(:,:,i),t);
    v(:,:,i) = spline(tk,vk(:,:,i),t);
    a(:,:,i) = spline(tk,ak(:,:,i),t); 
end


figure(1)
colors = get(gca,'colororder');
colors = [colors; [1,0,0];[0,1,0];[0,0,1]];


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
        pause(0.2)
    end
    clf
    pause(0.1)
end


diff = p(:,:,1) - p(:,:,2);
dist = sqrt(sum(diff.^2,1));

figure(5)
plot(t,dist, 'LineWidth',1.5);
grid on;
xlabel('t [s]')
ylabel('distance [m]');