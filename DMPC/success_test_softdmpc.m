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
k_hor = 15;
success = 1;
N_vector = 10:5:50; % number of vehicles
trials = 100;

% Variables for ellipsoid constraint
order = 4; % choose between 2 or 4 for the order of the super ellipsoid
rmin = 0.5; % X-Y protection radius for collisions
c = 1.5; % make this one for spherical constraint
E = diag([1,1,c]);
E1 = E^(-1);
E2 = E^(-order);

% Workspace boundaries
pmin = [-2.5,-2.5,0.2];
pmax = [2.5,2.5,2.2];

% Minimum distance between vehicles in m
rmin_init = 0.91;

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
tol = 2;
fail = 0;
error_tol = 0.05; % 5cm destination tolerance

Delta = getDeltaMat(k_hor); 

for k = 1:k_hor
    A_init = Aux*A_init;
    A_initp = [A_initp; A_init(1:3,:)];  
end

% Start Test

for q = 1:length(N_vector)
    N = N_vector(q);
    for r = 1:trials
        fprintf("Doing trial #%i with %i vehicles\n",r,N)
        % Initial positions
        [po,pf] = randomTest(N,pmin,pmax,rmin_init);
        
        %DMPC
        l = [];
        success = 0; %check if QP was feasible
        at_goal = 0; %At the end of solving, makes sure every agent arrives at the goal
        error_tol = 0.05; % 5cm destination tolerance
        violation(q,r) = 0; % checks if violations occured at end of algorithm

        % Penalty matrices when there're predicted collisions
        Q = 1000;
        S = 100;
        tries(q,r) = 1;
        failed_goal(q,r) = 0;
        t_start = tic;
        while tries(q,r) <= 1 && ~at_goal
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
                        [pi,vi,ai,success] = solveSoftDMPC(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta,Q,S,E1,E2,order); 
                    end
                    if ~success
                        break;
                    end
                    new_l(:,:,n) = pi;
                    pk(:,k,n) = pi(:,1);
                    vk(:,k,n) = vi(:,1);
                    ak(:,k,n) = ai(:,1);
                end
                if ~success
                    tries(q,r) = tries(q,r) + 1;
                    Q = Q+50;
                    break;
                end
                l = new_l;
            end
            if ~success
                continue
            end
            pass = ReachedGoal(pk,pf,K,error_tol);
            if success && pass
                at_goal = 1;
            elseif success && ~pass
                failed_goal(q,r) = failed_goal(q,r) + 1;
                tries(q,r) = tries(q,r) + 1;
                Q = Q+100;
            end
        end

        if success && at_goal      
            % Check if collision constraints were not violated
            for i = 1:N
                for j = 1:N
                    if(i~=j)
                        differ = E1*(pk(:,:,i) - pk(:,:,j));
                        dist = (sum(differ.^order,1)).^(1/order);
                        if min(dist) < (rmin - 0.05)
                            [value,index] = min(dist);
                            violation(q,r) = 1;
                        end
                    end
                end
            end
            
            for i = 1:N
                p(:,:,i) = spline(tk,pk(:,:,i),t);
                v(:,:,i) = spline(tk,vk(:,:,i),t);
                a(:,:,i) = spline(tk,ak(:,:,i),t); 
            end
            t_dmpc(q,r) = toc(t_start);
            totdist_dmpc(q,r) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
            
            for i = 1:N
                diff_goal = p(:,:,i) - repmat(pf(:,:,i),length(t),1)';
                dist_goal = sqrt(sum(diff_goal.^2,1));
                time_index(i) = find(dist_goal > 0.05,1,'last') + 1;
            end
            traj_time(q,r) = max(time_index)*Ts;
        else
            t_dmpc(q,r) = nan;
            totdist_dmpc(q,r) = nan;
%             save(['Fail_' num2str(fail)]);
            fail = fail + 1;
            traj_time(q,r) = nan;
        end
        success_dmpc(q,r) = success && at_goal && ~violation(q,r);
    end
end
fprintf("Finished! \n")
%% Post-Processing

% Probability of success plots
prob_dmpc = sum(success_dmpc,2)/trials;
figure(1)
grid on;
hold on;
ylim([0,1.05])
plot(N_vector,prob_dmpc,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Success Probability');

% Computation time
tmean_dmpc = nanmean(t_dmpc,2);
tstd_dmpc = nanstd(t_dmpc,1,2);
figure(2)
grid on;
hold on;
errorbar(N_vector,tmean_dmpc,tstd_dmpc,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Average Computation time [s]');

% Completion time
tmean_traj = nanmean(traj_time,2);
tstd_traj = nanstd(traj_time,1,2);
figure(3)
grid on;
hold on;
errorbar(N_vector,tmean_traj,tstd_traj,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Average Time for Transition [s]');

% Failure analysis
violation_num = sum(violation,2);
goal_num = sum(failed_goal,2);
total_num = sum(~success_dmpc,2);
infes_num = abs(total_num - goal_num);
figure(4)
grid on;
hold on;
bar(N_vector,[infes_num violation_num goal_num],'stacked');
xlabel('Number of Vehicles');
ylabel(['Number of failed trials (out of ' ,num2str(trials), ')']);
legend('Infeasibility','Collisions','Incomplete Trajectory')
