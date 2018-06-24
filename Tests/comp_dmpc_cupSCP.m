clc
clear all
close all
warning('off','all')

% Time settings and variables
T = 15; % Trajectory final time
h = 0.2; % time step duration
tk = 0:h:T;
K = T/h + 1; % number of time steps
Ts = 0.01; % period for interpolation @ 100Hz
t = 0:Ts:T; % interpolated time vector
k_hor = 15; % horizon length (currently set to 3s)
N_vector = 2:2:4; % number of vehicles
trials = 2; % number os trails per number of vehicles

% Variables for ellipsoid constraint
order = 2; % choose between 2 or 4 for the order of the super ellipsoid
rmin = 0.5; % X-Y protection radius for collisions
c = 1.5; % make this one for spherical constraint
E = diag([1,1,c]);
E1 = E^(-1);
E2 = E^(-order);

% Workspace boundaries
pmin = [-1.0,-1.0,0.2];
pmax = [1.0,1.0,2.2];

% Minimum distance between vehicles in m
rmin_init = 0.75;

% Maximum acceleration in m/s^2
alim = 0.5;

% Some pre computations DMPC
A = getPosMat(h,k_hor);
Aux = [1 0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];
A_initp = [];
A_init = eye(6);
fail = 0;
error_tol = 0.05; % 5cm destination tolerance

Delta = getDeltaMat(k_hor); 

for k = 1:k_hor
    A_init = Aux*A_init;
    A_initp = [A_initp; A_init(1:3,:)];  
end

b = [h^2/2*eye(3);
     h*eye(3)];

prev_row = zeros(6,3*K); % For the first iteration of constructing matrix Ain
A_p = zeros(3*(K-1),3*K);
A_v = zeros(3*(K-1),3*K);
idx=1;
% Build matrix to convert acceleration to position
for k = 1:(K-1)
    add_b = [zeros(size(b,1),size(b,2)*(k-1)) b zeros(size(b,1),size(b,2)*(K-k))];
    new_row = Aux*prev_row + add_b;   
    A_p(idx:idx+2,:) = new_row(1:3,:);
    A_v(idx:idx+2,:) = new_row(4:6,:);
    prev_row = new_row;
    idx = idx+3;
end

% Start Test

for q = 1:length(N_vector)
    N = N_vector(q);
    for r = 1:trials
        fprintf("Doing trial #%i with %i vehicles\n",r,N)
        % Initial positions
        [po,pf] = randomTest(N,pmin,pmax,rmin_init);
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
        %SoftDMPC with bound on relaxation variable
        
        % Variables for ellipsoid constraint
        order = 2; % choose between 2 or 4 for the order of the super ellipsoid
        rmin = 0.5; % X-Y protection radius for collisions
        c = 1.5; % make this one for spherical constraint
        E = diag([1,1,c]);
        E1 = E^(-1);
        E2 = E^(-order);
        term4 = -5*10^4;
        
        l = [];
        feasible(q,r) = 0; %check if QP was feasible
        error_tol = 0.05; % 5cm destination tolerance
        violation(q,r) = 0; % checks if violations occured at end of algorithm

        % Penalty matrices when there're predicted collisions
        Q = 1000;
        S = 100;
        failed_goal(q,r) = 0;
        outbound(q,r) = 0;
        t_start = tic;
        
        for k = 1:K
            for n = 1:N
                if k==1
                    poi = po(:,:,n);
                    pfi = pf(:,:,n);
                    [pi,vi,ai] = initDMPC(poi,pfi,h,k_hor,K);
                    feasible(q,r) = 1;
                else
                    pok = pk(:,k-1,n);
                    vok = vk(:,k-1,n);
                    aok = ak(:,k-1,n);
                    [pi,vi,ai,feasible(q,r),outbound(q,r)] = solveSoftDMPCbound(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta,Q,S,E1,E2,order,term4); 
                end
                if ~feasible(q,r)
                    break;
                end
                new_l(:,:,n) = pi;
                pk(:,k,n) = pi(:,1);
                vk(:,k,n) = vi(:,1);
                ak(:,k,n) = ai(:,1);
            end
            if ~feasible(q,r)
                break;
            end
            l = new_l;
        end
        if feasible(q,r)
            pass = ReachedGoal(pk,pf,K,error_tol,N);
            if  ~pass
                failed_goal(q,r) = failed_goal(q,r) + 1;
            end
        end

        if feasible(q,r) && ~failed_goal(q,r)       
            for i = 1:N
                p(:,:,i) = spline(tk,pk(:,:,i),t);
                v(:,:,i) = spline(tk,vk(:,:,i),t);
                a(:,:,i) = spline(tk,ak(:,:,i),t); 
            end
            % Check if collision constraints were not violated
            for i = 1:N
                for j = 1:N
                    if(i~=j)
                        differ = E1*(p(:,:,i) - p(:,:,j));
                        dist = (sum(differ.^order,1)).^(1/order);
                        if min(dist) < (rmin - 0.05)
                            [value,index] = min(dist);
                            violation(q,r) = 1;
                        end
                    end
                end
            end
            t_dmpc(q,r) = toc(t_start);
            totdist_dmpc(q,r) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
            
            for i = 1:N
                diff_goal = p(:,:,i) - repmat(pf(:,:,i),length(t),1)';
                dist_goal = sqrt(sum(diff_goal.^2,1));
                hola = find(dist_goal >= 0.05,1,'last');
                if isempty(hola)
                    time_index(i) = 0;
                else
                    time_index(i) = hola + 1;
                end
            end
            traj_time(q,r) = max(time_index)*Ts;
        else
            t_dmpc(q,r) = nan;
            totdist_dmpc(q,r) = nan;
            traj_time(q,r) = nan;
        end
        success_dmpc(q,r) = feasible(q,r) && ~failed_goal(q,r) && ~violation(q,r);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cup-SCP
        p = [];
        v = [];
        a = [];
        t_start = tic;
        % Solve SCP
        [pk,vk,ak,success_cup(q,r)] = solveCupSCP(po,pf,h,K,N,pmin,pmax,rmin,alim,A_p,A_v,E1,E2,order);
        
        if (success_cup(q,r))
            % Interpolate solution with a 100Hz sampling
            for i = 1:N
                p(:,:,i) = spline(tk,pk(:,:,i),t);
                v(:,:,i) = spline(tk,vk(:,:,i),t);
                a(:,:,i) = spline(tk,ak(:,:,i),t); 
            end
            t_cup(q,r) = toc(t_start);
            totdist_cup(q,r) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));        
        else
            t_cup(q,r) = nan;
            totdist_cup(q,r) = nan;
        end    
    end
end
fprintf("Finished! \n")
save('comp_dmpc_cupSCP3')
%% Post-Processing

% Probability of success plots
prob_dmpc = sum(success_dmpc,2)/trials;
prob_cup = sum(success_cup,2)/trials;
figure(1)
grid on;
hold on;
ylim([0,1.05])
plot(N_vector,prob_dmpc,'Linewidth',2);
plot(N_vector,prob_cup,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Success Probability');
legend('DMPC','cup-SCP');

% Computation time
tmean_dmpc = nanmean(t_dmpc,2);
tstd_dmpc = nanstd(t_dmpc,1,2);
tmean_cup = nanmean(t_cup,2);
tstd_cup = nanstd(t_cup,1,2);
figure(2)
grid on;
hold on;
errorbar(N_vector,tmean_dmpc,tstd_dmpc,'Linewidth',2);
errorbar(N_vector,tmean_cup,tstd_cup,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Average Computation time [s]');
legend('DMPC','cup-SCP');

% Percentage increase/decrease on travelled dist of dmpc wrt dec
% Positive number means that dmpc path was longer
diff_dist = (totdist_dmpc-totdist_cup)./totdist_dmpc;
avg_diff = nanmean(diff_dist,2);
std_diff = nanstd(diff_dist,1,2);
figure(3)
errorbar(N_vector,100*avg_diff,100*std_diff,'Linewidth',2);
grid on;
xlabel('Number of Vehicles');
ylabel('Average % increase');
title('Percentage increase on travelled distance of DMPC wrt cup-SCP');


% Failure analysis
violation_num = sum(violation,2);
goal_num = sum(failed_goal,2);
infes_num = sum(~feasible,2);
figure(4)
grid on;
hold on;
bar(N_vector,[infes_num violation_num goal_num],'stacked');
xlabel('Number of Vehicles');
ylabel(['Number of failed trials (out of ' ,num2str(trials), ')']);
legend('Infeasibility','Collisions','Incomplete Trajectory')
