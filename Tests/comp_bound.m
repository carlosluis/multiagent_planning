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
k_hor = 15; % horizon length (currently set to 3s)
N_vector = 10:10:80; % number of vehicles
trials = 50; % number os trails per number of vehicles
fail2 = 0;
fail4 = 0;

% Workspace boundaries
pmin = [-2.5,-2.5,0.2];
pmin_gen = [-2.4,-2.4,0.3];
pmax = [2.5,2.5,2.2];
pmax_gen = [2.4,2.4,2.1];

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
        
        %SoftDMPC with second-order ellipsoid constraint
        
        % Variables for ellipsoid constraint
        order = 2; % choose between 2 or 4 for the order of the super ellipsoid
        rmin = 0.5; % X-Y protection radius for collisions
        c = 1.5; % make this one for spherical constraint
        E = diag([1,1,c]);
        E1 = E^(-1);
        E2 = E^(-order);
        term2 = -1*10^4;
        
        l = [];
        feasible2(q,r) = 0; %check if QP was feasible
        error_tol = 0.05; % 5cm destination tolerance
        violation2(q,r) = 0; % checks if violations occured at end of algorithm

        % Penalty matrices when there're predicted collisions
        Q = 1000;
        S = 100;
        failed_goal2(q,r) = 0;
        outbound2(q,r) = 0;
        t_start = tic;
        
        for k = 1:K
            for n = 1:N
                if k==1
                    poi = po(:,:,n);
                    pfi = pf(:,:,n);
                    [pi,vi,ai] = initDMPC(poi,pfi,h,k_hor,K);
                    feasible2(q,r) = 1;
                else
                    pok = pk(:,k-1,n);
                    vok = vk(:,k-1,n);
                    aok = ak(:,k-1,n);
                    [pi,vi,ai,feasible2(q,r),outbound2(q,r)] = solveSoftDMPCrepair(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta,Q,S,E1,E2,order,term2); 
                end
                if ~feasible2(q,r)
                    save(['Fail2_' num2str(fail2)]);
                    fail2 = fail2 + 1;
                    break;
                end
                new_l(:,:,n) = pi;
                pk(:,k,n) = pi(:,1);
                vk(:,k,n) = vi(:,1);
                ak(:,k,n) = ai(:,1);
            end
            if ~feasible2(q,r)
                break;
            end
            l = new_l;
        end
        if feasible2(q,r)
            pass = ReachedGoal(pk,pf,K,error_tol,N);
            if  ~pass
                failed_goal2(q,r) = failed_goal2(q,r) + 1;
            end
        end

        if feasible2(q,r) && ~failed_goal2(q,r)               
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
                            violation2(q,r) = 1;
                        end
                    end
                end
            end
            t_dmpc2(q,r) = toc(t_start);
            totdist_dmpc2(q,r) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
            
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
            traj_time2(q,r) = max(time_index)*Ts;
        else
            t_dmpc2(q,r) = nan;
            totdist_dmpc2(q,r) = nan;
            traj_time2(q,r) = nan;
        end
        success_dmpc2(q,r) = feasible2(q,r) && ~failed_goal2(q,r) && ~violation2(q,r);
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
        %SoftDMPC with bound on relaxation variable
        
        % Variables for ellipsoid constraint
        order = 2; % choose between 2 or 4 for the order of the super ellipsoid
        rmin = 0.5; % X-Y protection radius for collisions
        c = 1.5; % make this one for spherical constraint
        E = diag([1,1,c]);
        E1 = E^(-1);
        E2 = E^(-order);
        term4 = -1*10^4;
        
        l = [];
        feasible4(q,r) = 0; %check if QP was feasible
        error_tol = 0.05; % 5cm destination tolerance
        violation4(q,r) = 0; % checks if violations occured at end of algorithm

        % Penalty matrices when there're predicted collisions
        Q = 1000;
        S = 100;
        failed_goal4(q,r) = 0;
        outbound4(q,r) = 0;
        t_start = tic;
        
        for k = 1:K
            for n = 1:N
                if k==1
                    poi = po(:,:,n);
                    pfi = pf(:,:,n);
                    [pi,vi,ai] = initDMPC(poi,pfi,h,k_hor,K);
                    feasible4(q,r) = 1;
                else
                    pok = pk(:,k-1,n);
                    vok = vk(:,k-1,n);
                    aok = ak(:,k-1,n);
                    [pi,vi,ai,feasible4(q,r),outbound4(q,r)] = solveSoftDMPCbound(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta,Q,S,E1,E2,order,term4); 
                end
                if ~feasible4(q,r)
                    break;
                end
                new_l(:,:,n) = pi;
                pk(:,k,n) = pi(:,1);
                vk(:,k,n) = vi(:,1);
                ak(:,k,n) = ai(:,1);
            end
            if ~feasible4(q,r)
                save(['Fail4_' num2str(fail4)]);
                fail4 = fail4 + 1;
                break;
            end
            l = new_l;
        end
        if feasible4(q,r)
            pass = ReachedGoal(pk,pf,K,error_tol,N);
            if  ~pass
                failed_goal4(q,r) = failed_goal4(q,r) + 1;
            end
        end

        if feasible4(q,r) && ~failed_goal4(q,r)       
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
                            violation4(q,r) = 1;
                        end
                    end
                end
            end
            t_dmpc4(q,r) = toc(t_start);
            totdist_dmpc4(q,r) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
            
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
            traj_time4(q,r) = max(time_index)*Ts;
        else
            t_dmpc4(q,r) = nan;
            totdist_dmpc4(q,r) = nan;
            traj_time4(q,r) = nan;
        end
        success_dmpc4(q,r) = feasible4(q,r) && ~failed_goal4(q,r) && ~violation4(q,r);
    end
end
fprintf("Finished! \n")
save('comp_bound_1')
%% Post-Processing
close all

% Volumen of arena
V = (pmax(1)-pmin(1))*(pmax(2)-pmin(2))*(pmax(3)-pmin(3));

% Probability of success plots
prob_dmpc2 = sum(success_dmpc2,2)/trials;
prob_dmpc4 = sum(success_dmpc4,2)/trials;

figure(1)
grid on;
hold on;
ylim([0,1.05])
plot(N_vector/V,prob_dmpc2,'Linewidth',2);
plot(N_vector/V,prob_dmpc4,'Linewidth',2);
xlabel('Workspace density [agent/m^3]');
ylabel('Success Probability');
legend('Repair','Bound')

% Computation time
tmean_dmpc2 = nanmean(t_dmpc2,2);
tstd_dmpc2 = nanstd(t_dmpc2,1,2);
tmean_dmpc4 = nanmean(t_dmpc4,2);
tstd_dmpc4 = nanstd(t_dmpc4,1,2);
figure(2)
grid on;
hold on;
errorbar(N_vector/V,tmean_dmpc2,tstd_dmpc2,'Linewidth',2);
errorbar(N_vector/V,tmean_dmpc4,tstd_dmpc4,'Linewidth',2);
xlabel('Workspace density [agent/m^3]');
ylabel('Average Computation time [s]');
legend('Repair','Bound')

% Completion time
tmean_traj2 = nanmean(traj_time2,2);
tstd_traj2 = nanstd(traj_time2,1,2);
tmean_traj4 = nanmean(traj_time4,2);
tstd_traj4 = nanstd(traj_time4,1,2);
figure(3)
grid on;
hold on;
errorbar(N_vector/V,tmean_traj2,tstd_traj2,'Linewidth',2);
errorbar(N_vector/V,tmean_traj4,tstd_traj4,'Linewidth',2);
xlabel('Workspace density [agent/m^3]');
ylabel('Average Time for Transition [s]');
legend('Repair','Bound')

% Failure analysis
violation_num2 = sum(violation2,2);
goal_num2 = sum(failed_goal2,2);
infes_num2 = sum(~feasible2,2);
total_num2 = sum(violation_num2) + sum(goal_num2) + sum(infes_num2);

violation_num4 = sum(violation4,2);
goal_num4 = sum(failed_goal4,2);
infes_num4 = sum(~feasible4,2);
total_num4 = sum(violation_num4) + sum(goal_num4) + sum(infes_num4);

StackData2 = [infes_num2/trials*100 violation_num2/trials*100 goal_num2/trials*100];
StackData4 = [infes_num4/trials*100 violation_num4/trials*100 goal_num4/trials*100];

figure(4)
h2=bar(N_vector/V-0.01,StackData2,'stacked','BarWidth',0.3);
grid on;
hold on;
h4 = bar(N_vector/V+0.01,StackData4,'stacked','BarWidth',0.3);
myC2= summer(size(StackData4,2));
myC4= winter(size(StackData4,2));
for k = 1:size(StackData4,2)
    set(h2(k),'facecolor',myC2(k,:));
    set(h4(k),'facecolor',myC4(k,:));
end
xticks(N_vector/V);
% xlim([min(N_vector/V)-1, max(N_vector/V)+1])
xlabel('Workspace density [agent/m^3]');
ylabel(['Failure rate % (out of ' ,num2str(trials), ')']);
legend('Infeasibility Rep','Collisions Rep','Incomplete Trajectory Rep',...
       'Infeasibility Bnd','Collisions Bnd','Incomplete Trajectory Bnd')
   
% Print various statistics
fprintf("--------------------------- \nREPAIR TEST STATISTICS \n");
fprintf("--------------------------- \n")
fprintf("Probability of success across all tests --> %.2f%% (%d out of %d tests) \n",...
    sum(prob_dmpc2)/length(prob_dmpc2)*100, sum(sum(success_dmpc2)), length(N_vector)*trials)
fprintf("Percentage failure due to infeasibility --> %.2f%% (%d out of %d failures) \n",...
    sum(infes_num2)/total_num2*100, sum(infes_num2), total_num2)

fprintf("Percentage failure due to collisions --> %.2f%% (%d out of %d failures) \n",...
    sum(violation_num2)/total_num2*100, sum(violation_num2), total_num2)

fprintf("Percentage failure due to not reaching goal --> %.2f%% (%d out of %d failures) \n",...
    sum(goal_num2)/total_num2*100, sum(goal_num2), total_num2)

fprintf("--------------------------- \nBOUND TEST STATISTICS \n");
fprintf("--------------------------- \n")
fprintf("Probability of success across all tests --> %.2f%% (%d out of %d tests) \n",...
    sum(prob_dmpc4)/length(prob_dmpc4)*100, sum(sum(success_dmpc4)), length(N_vector)*trials)
fprintf("Percentage failure due to infeasibility --> %.2f%% (%d out of %d failures) \n",...
    sum(infes_num4)/total_num4*100, sum(infes_num4), total_num4)

fprintf("Percentage failure due to collisions --> %.2f%% (%d out of %d failures) \n",...
    sum(violation_num4)/total_num4*100, sum(violation_num4), total_num4)

fprintf("Percentage failure due to not reaching goal --> %.2f%% (%d out of %d failures) \n",...
    sum(goal_num4)/total_num4*100, sum(goal_num4), total_num4)


