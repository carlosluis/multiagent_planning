clc
clear all
close all
warning('off','all')

% Time settings and variables
T = 15; % Trajectory final time
h = 0.2; % time step duration
tk = 0:h:T;
K = T/h + 1; % number of time steps
k_hor = 15; % horizon length (currently set to 3s)
Ts = 0.01; % period for interpolation @ 100Hz
t = 0:Ts:T; % interpolated time vector
N_vector = 2:2:20; % number of vehicles
trials = 50; % number os trails per number of vehicles

% Variables for ellipsoid constraint
order = 2; % choose between 2 or 4 for the order of the super ellipsoid
rmin = 0.35; % X-Y protection radius for collisions
c = 2.0; % make this one for spherical constraint
E = diag([1,1,c]);
E1 = E^(-1);
E2 = E^(-order);

% Workspace boundaries
pmin = [-1.0,-1.0,0.2];
pmax = [1.0,1.0,2.2];

% Minimum distance between vehicles in m
rmin_init = 0.75;

% Maximum acceleration in m/s^2
alim = 1.0;

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
fail = 0;

for q = 1:length(N_vector)
    N = N_vector(q);
    for r = 1:trials
        fprintf("Doing trial #%i with %i vehicles\n",r,N)
        % Initial positions
        [po,pf] = randomTest(N,pmin,pmax,rmin_init);
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
        %SoftDMPC with bound on relaxation variable
        T = 0;
        l = [];
        p = [];
        v = [];
        a = [];
        pk = [];
        vk = [];
        ak = [];
        pred = [];
        coll(q,r) = 0;
        term = -5*10^4;
    
        feasible(q,r) = 0; %check if QP was feasible
        error_tol = 0.01; % 5cm destination tolerance
        violation(q,r) = 0; % checks if violations occured at end of algorithm

        % Penalty matrices when there're predicted collisions
        Q = 1000;
        S = 100;
        failed_goal(q,r) = 0;
        outbound(q,r) = 0;
        reached_goal = 0;
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
                    [pi,vi,ai,feasible(q,r),outbound(q,r),coll(q,r)] = solveSoftDMPCbound2(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta,Q,S,E1,E2,order,term); 
                end
                if (~feasible(q,r) || outbound(q,r) || coll(q,r)) %problem was infeasible, exit and retry
                    break;
                end
                new_l(:,:,n) = pi;
                pk(:,k,n) = pi(:,1);
                vk(:,k,n) = vi(:,1);
                ak(:,k,n) = ai(:,1);
            end
            if ~feasible(q,r)
                save(['Fail_' num2str(fail)]); 
                fail = fail + 1;
                break;
            end
            l = new_l;
            pred(:,:,:,k) = l;
        end
        
        if feasible(q,r)
            reached_goal = ReachedGoal(pk,pf,k,error_tol,N);       
            if ~reached_goal
                failed_goal(q,r) = 1;
            end   
        end
        
        if feasible(q,r) && ~failed_goal(q,r)             
            % Interpolate @ 100 Hz
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
                hola = find(dist_goal >= error_tol,1,'last');
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
        pk = [];
        vk = [];
        ak = [];
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dec-iSCP
        % Empty list of obstacles
        l = [];
        pk = [];
        vk = [];
        ak = [];
        p = [];
        v = [];
        a = [];
        
        % DEC-ISCP
        t_start = tic; 
        for i = 1:N 
            poi = po(:,:,i);
            pfi = pf(:,:,i);
            [pi, vi, ai,success] = singleiSCP(poi,pfi,h,K,pmin,pmax,rmin,alim,l,A_p,A_v,E1,E2,order);
            if ~success
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
            t_dec(q,r) = toc(t_start);
            totdist_dec(q,r) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
        
        else
            t_dec(q,r) = nan;
            totdist_dec(q,r) = nan;
        end
        success_dec(q,r) = success;
        
        
    end
end
fprintf("Finished! \n")
save('comp_all_fixedT_1')
%% Post-Processing

% Probability of success plots
prob_dmpc = sum(success_dmpc,2)/trials;
prob_cup = sum(success_cup,2)/trials;
prob_dec = sum(success_dec,2)/trials;
figure(1)
grid on;
hold on;
ylim([0,1.05])
plot(N_vector,prob_cup,'Linewidth',2);
plot(N_vector,prob_dec','Linewidth',2);
plot(N_vector,prob_dmpc,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Success Probability');
legend('cup-SCP','dec-iSCP','DMPC');

% Computation time
tmean_dmpc = nanmean(t_dmpc,2);
tstd_dmpc = nanstd(t_dmpc,1,2);
tmean_cup = nanmean(t_cup,2);
tstd_cup = nanstd(t_cup,1,2);
tmean_dec = nanmean(t_dec,2);
tstd_dec = nanstd(t_dec,1,2);
figure(2)
grid on;
hold on;
plot(N_vector, tmean_cup,'LineWidth',2);
plot(N_vector, tmean_dec,'LineWidth',2);
plot(N_vector, tmean_dmpc,'LineWidth',2);
% errorbar(N_vector,tmean_cup,tstd_cup,'Linewidth',2);
% errorbar(N_vector,tmean_dmpc,tstd_dmpc,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Average Computation time [s]');
legend('cup-SCP','dec-iSCP','DMPC');

% Average travelled distance
avg_dist_dmpc = nanmean(totdist_dmpc,2);
avg_dist_cup = nanmean(totdist_cup,2);
avg_dist_dec = nanmean(totdist_dec,2);
figure(3)
plot(N_vector, avg_dist_cup,'LineWidth', 3);
hold on;
grid on;
plot(N_vector, avg_dist_dec,'LineWidth', 3);
plot(N_vector, avg_dist_dmpc,'LineWidth', 3);
xlabel('Number of Vehicles');
ylabel('Total Travelled Distance [m]');
legend('cup-SQP','dec-iSCP','DMPC');

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
