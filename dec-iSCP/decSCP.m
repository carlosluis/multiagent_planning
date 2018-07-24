clc
clear all
close all

% Time settings and variables
T = 7; % Trajectory final time
h = 0.2; % time step duration
tk = 0:h:T;
K = T/h + 1; % number of time steps
Ts = 0.01; % period for interpolation @ 100Hz
t = 0:Ts:T; % interpolated time vector
success = 1;
N = 2; % number of vehicles

% Variables for ellipsoid constraint
order = 2; % choose between 2 or 4 for the order of the super ellipsoid
rmin = 0.31; % X-Y protection radius for collisions
c = 2.0; % make this one for spherical constraint
E = diag([1,1,c]);
E1 = E^(-1);
E2 = E^(-order);

% Workspace boundaries
pmin = [-1.0,-1.0,0.2];
pmax = [1.0,1.0,2.2];

rmin_init = 0.75;

% Initial positions
% [po,pf] = randomTest(N,pmin,pmax,rmin_init);

po1 = [-1.0, -1.0, 1.0];
po2 = [1.0, 1.0, 1.0];

po = cat(3,po1,po2);
pf = cat(3,po2,po1);

%% Some Precomputations
l = [];
p = [];
v = [];
a = [];
pk = [];
vk = [];
ak = [];
reached_goal = 0;
violation = 0;
error_tol = 0.01;
% Kinematic model A,b matrices
A = [1 0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

b = [h^2/2*eye(3);
     h*eye(3)];
 
prev_row = zeros(6,3*K); % For the first iteration of constructing matrix Ain
A_p = zeros(3*(K-1),3*K);
A_v = zeros(3*(K-1),3*K);
idx=1;
% Build matrix to convert acceleration to position
for k = 1:(K-1)
    add_b = [zeros(size(b,1),size(b,2)*(k-1)) b zeros(size(b,1),size(b,2)*(K-k))];
    new_row = A*prev_row + add_b;   
    A_p(idx:idx+2,:) = new_row(1:3,:);
    A_v(idx:idx+2,:) = new_row(4:6,:);
    prev_row = new_row;
    idx = idx+3;
end

% Empty list of obstacles
l = [];

% Maximum acceleration in m/s^2
alim = 1.0;

tic %measure the time it gets to solve the optimization problem
for i = 1:N 
    poi = po(:,:,i);
    pfi = pf(:,:,i);
    [pi, vi, ai,success] = singleiSCP(poi,pfi,h,K,pmin,pmax,rmin,alim,l,A_p,A_v,E1,E2,order);
    if ~success
        fprintf('Failed solving for vehicle %i\n',i);
        break;
    end
    l = cat(3,l,pi);
    pk(:,:,i) = pi;
    vk(:,:,i) = vi;
    ak(:,:,i) = ai;
    
    % Interpolate solution with a 100Hz sampling
    p(:,:,i) = spline(tk,pi,t);
    v(:,:,i) = spline(tk,vi,t);
    a(:,:,i) = spline(tk,ai,t);
end

if success
    reached_goal = ReachedGoal(pk,pf,K,error_tol,N);
    % Check if collision constraints were not violated
    for i = 1:N
        for j = 1:N
            if(i~=j)
                differ = E1*(p(:,:,i) - p(:,:,j));
                dist = (sum(differ.^order,1)).^(1/order);
                if min(dist) < (rmin-0.01)
                    [value,index] = min(dist);
                    violation = 1;
                    fprintf("Collision constraint violated by %.2fcm: vehicles %i and %i @ k = %i \n", (rmin -value)*100,i,j,index)
                end
            end
        end
    end
end

pass = success && reached_goal && ~violation
toc

%%
L = length(t);
colors = distinguishable_colors(N);
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
colors = distinguishable_colors(N);
for i = 1:N
    figure(1);
    h_plot(i) = plot3(p(1,:,i), p(2,:,i), p(3,:,i), 'LineWidth',1.5,...
                'Color',colors(i,:));
    h_label{i} = ['Vehicle #' num2str(i)];
    hold on;
    grid on;
    xlim([-4,4])
    ylim([-4,4])
    zlim([0,3.5])
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
plot(t,rmin*ones(length(tk),1),'--r','LineWidth',1.5);
legend(h_plot,h_label);