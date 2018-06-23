clc
clear all
close all
warning('off','all')

% Time settings and variables
T = 20; % Trajectory final time
h = 0.2; % time step duration
tk = 0:h:T;
K = T/h + 1; % number of time steps
Ts = 0.01; % period for interpolation @ 100Hz
t = 0:Ts:T; % interpolated time vector
k_hor = 15; % horizon length

% Variables for ellipsoid constraint
order = 2; % choose between 2 or 4 for the order of the super ellipsoid
rmin = 0.5; % X-Y protection radius for collisions
c = 1.5; % make this one for spherical constraint
E = diag([1,1,c]);
E1 = E^(-1);
E2 = E^(-order);

N = 80; % number of vehicles

% Workspace boundaries
pmin = [-2.5,-2.5,0.2];
pmax = [2.5,2.5,2.2];

% % Workspace boundaries
% pmin = [-5,-5,0.2];
% pmax = [5,5,5];

% Minimum distance between vehicles in m
rmin_init = 0.75;

% Initial positions
[po,pf] = randomTest(N,pmin,pmax,rmin_init);

% % Initial positions
% po1 = [1.501,1.5,1.5];
% po2 = [-1.5,-1.5,1.5];
% po3 = [-1.5,1.5,1.5];
% po4 = [1.5,-1.5,1.5];
% po = cat(3,po1,po2,po3,po4);
% 
% % Final positions
% pf1 = [-1.5,-1.5,1.5];
% pf2 = [1.5,1.5,1.5];
% pf3 = [1.5,-1.5,1.5];
% pf4 = [-1.5,1.5,1.5];
% pf  = cat(3,pf1,pf2,pf3,pf4);

%% Solving the problem
l = [];
p = [];
v = [];
a = [];
pk = [];
vk = [];
ak = [];
success = 0; %check if QP was feasible
at_goal = 0; %At the end of solving, makes sure every agent arrives at the goal
error_tol = 0.05; % 5cm destination tolerance
violation = 0; % checks if violations occured at end of algorithm
outbound = 0;
term = -5*10^4;

% Penalty matrices when there're predicted collisions
Q = 1000;
S = 100;

% Maximum acceleration in m/s^2
alim = 0.5;

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

failed_goal = 0; %how many times the algorithm failed to reach goal
tries = 1; % how many iterations it took the DMPC to find a solution
tic % measure the time it gets to solve the optimization problem
pred = [];
moreconstr = [];
for k = 1:K
        for n = 1:N
        if k==1
            poi = po(:,:,n);
            pfi = pf(:,:,n);
            [pi,vi,ai] = initDMPC(poi,pfi,h,k_hor,K);
            success = 1;
        else
            pok = pk(:,k-1,n);
            vok = vk(:,k-1,n);
            aok = ak(:,k-1,n);
            [pi,vi,ai,success,outbound] = solveSoftDMPCbound(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta,Q,S,E1,E2,order,term); 
        end
        if (~success || outbound) %problem was infeasible, exit and retry
            break;
        end
        new_l(:,:,n) = pi;
        pk(:,k,n) = pi(:,1);
        vk(:,k,n) = vi(:,1);
        ak(:,k,n) = ai(:,1);
    end
    if ~success %Heuristic: increase Q, make init more slowly,
        if outbound
            fprintf("Failed - problem unfeasible, vehicle couldn't stay in workspace @ k_T = %i, n = %i: trial #%i\n",k,n,tries-1)
        else
            fprintf("Failed - problem unfeasible @ k_T = %i, n = %i: trial #%i\n",k,n,tries-1)
        end
        break;
    end
    l = new_l;
    pred(:,:,:,k) = l;
end
pass = 0;
if success
pass = ReachedGoal(pk,pf,K,error_tol,N); %check if agents reached goal
end
if success && pass
    at_goal = 1;
elseif success && ~pass %if not at goal, retry with more aggressive behaviour
    failed_goal = failed_goal + 1;
    fprintf("Failed - Did not reach goal: trial #%i\n",tries-1)
end
passed = success && at_goal %DMPC was successful or not      
toc

if passed
    
    % Interpolate for better resolution
    for i = 1:N
        p(:,:,i) = spline(tk,pk(:,:,i),t);
        v(:,:,i) = spline(tk,vk(:,:,i),t);
        a(:,:,i) = spline(tk,ak(:,:,i),t); 
    end
    
    % Check if collision constraints were not violated
    for i = 1:N
        for j = 1:N
            if(i~=j)
                differ = E1*(pk(:,:,i) - pk(:,:,j));
                dist = (sum(differ.^order,1)).^(1/order);
                if min(dist) < (rmin - 0.05)
                    [value,index] = min(dist);
                    violation = 1;
                    fprintf("Collision constraint violated by %.2fcm: vehicles %i and %i @ k = %i \n", (rmin -value)*100,i,j,index)
                end
            end
        end
    end

    if ~violation
        fprintf("No collisions found! Successful computation\n")
    end
    
    % Calculate how much time is required to complete the transition
    % within a 5cm margin of the goal

    for i = 1:N
        differ = p(:,:,i) - repmat(pf(:,:,i),length(t),1)';
        dist = sqrt(sum(differ.^2,1));
        hola = find(dist >= 0.05,1,'last');
        if isempty(hola)
            time_index(i) = 0;
        else
            time_index(i) = hola + 1;
        end
    end
    max_time_index = max(time_index);
    fprintf("The trajectory can be completed in %.2f seconds\n",max_time_index*Ts);
    totdist_dmpc = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
    fprintf("The sum of trajectory length is %.2f\n",totdist_dmpc);
end

%%
figure(1)
colors = distinguishable_colors(N);

set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'currentchar',' ')
while get(gcf,'currentchar')==' '
   
    for i = 1:N
    h_line(i) = animatedline('LineWidth',2,'Color',colors(i,:),'LineStyle',':');
    end
    for k = 1:K
        for i = 1:N
            clearpoints(h_line(i));
            addpoints(h_line(i),pred(1,:,i,k),pred(2,:,i,k),pred(3,:,i,k));     
            hold on;
            grid on;
            xlim([pmin(1),pmax(1)])
            ylim([pmin(2),pmax(2)])
            zlim([0,pmax(3)])
            plot3(pk(1,k,i),pk(2,k,i),pk(3,k,i),'o',...
                'LineWidth',2,'Color',colors(i,:));
            plot3(po(1,1,i), po(1,2,i), po(1,3,i),'^',...
                  'LineWidth',2,'Color',colors(i,:));
            plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i),'x',...
                  'LineWidth',2,'Color',colors(i,:));    
        end
    drawnow
    end
    clf
    pause(0.5)
end

%% Plotting
L = length(t);
colors = distinguishable_colors(N);
       
for i = 1:N
    figure(1);
    h_plot(i) = plot3(p(1,:,i), p(2,:,i), p(3,:,i), 'LineWidth',1.5,...
                'Color',colors(i,:));
    h_label{i} = ['Vehicle #' num2str(i)];
    hold on;
    grid on;
    xlim([pmin(1),pmax(1)])
    ylim([pmin(2),pmax(2)])
    zlim([0,pmax(3)])
    xlabel('x[m]')
    ylabel('y[m]');
    zlabel('z[m]')
    plot3(po(1,1,i), po(1,2,i), po(1,3,i),'x',...
                  'LineWidth',3,'Color',colors(i,:));
%     plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i),'x',...
%                   'LineWidth',5,'Color',colors(i,:)); 
    
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
%%
figure(6)
for i = 1:N
    for j = 1:N
        if(i~=j)
            differ = E1*(pk(:,:,i) - pk(:,:,j));
            dist = (sum(differ.^order,1)).^(1/order);
            plot(tk, dist, 'LineWidth',1.5);
            grid on;
            hold on;
            xlabel('t [s]')
            ylabel('Inter-agent distance [m]');
        end
    end
end
plot(tk,rmin*ones(length(tk),1),'--r','LineWidth',1.5);
% legend(h_plot,h_label);