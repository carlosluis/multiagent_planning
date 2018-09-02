clc
clear all
close all
warning('off','all')

% Time settings and variables
max_T = 30; % Trajectory final time
h = 0.2; % time step duration
max_K = max_T/h + 1; % number of time steps
k_hor = 15; % horizon length (currently set to 3s)
N_vector = 20:20:200; % number of vehicles
trials = 50; % number os trails per number of vehicles

% Variables for ellipsoid constraint
order = 2; % choose between 2 or 4 for the order of the super ellipsoid
c = 2.0; % make this one for spherical constraint
E = diag([1,1,c]);
E1 = E^(-1);
E2 = E^(-order);
fail = 0;

% Minimum distance between vehicles in m
rmin_init = 0.35;

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
 
b = [h^2/2*eye(3);
     h*eye(3)];
 
prev_row = zeros(6,3*k_hor); % For the first iteration of constructing matrix Ain
A_p_dmpc = [];
A_v_dmpc = [];
A_initp = [];
A_init = eye(6);
for k = 1:k_hor
    add_b = [zeros(size(b,1),size(b,2)*(k-1)) b zeros(size(b,1),size(b,2)*(k_hor-k))];
    new_row = Aux*prev_row + add_b;   
    A_p_dmpc = [A_p_dmpc; new_row(1:3,:)];
    A_v_dmpc = [A_v_dmpc; new_row(4:6,:)];
    prev_row = new_row;
    A_init = Aux*A_init;
    A_initp = [A_initp; A_init(1:3,:)];  
end

Delta = getDeltaMat(k_hor);

% Start Test
% We want to keep the density the same for all scenarios, 
% so we make it adaptive based on N

for q = 1:length(N_vector)
    N = N_vector(q);
    pmin = [-N^(1/3)/2 -N^(1/3)/2 0.2];
    pmax = [N^(1/3)/2 N^(1/3)/2 N^(1/3)+0.2];
    for r = 1:trials
        fprintf("Doing trial #%i with %i vehicles\n",r,N)
        % Initial positions
        [po,pf] = randomTest(N,pmin,pmax,rmin_init,E1,order);
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
        %SoftDMPC enforcing constraint at kctr = k
        T = 0;
        l = [];
        p = [];
        v = [];
        a = [];
        pk = [];
        vk = [];
        ak = [];
        coll(q,r) = 0;
        term = -5*10^4;
        
        % Variables for ellipsoid constraint
        rmin = 0.35; % X-Y protection radius for collisions
    
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
        k = 1;
        while ~reached_goal && k < max_K
            for n = 1:N
                if k==1
                    poi = po(:,:,n);
                    pfi = pf(:,:,n);
                    [pi,vi,ai] = initDMPC(poi,pfi,h,k_hor,max_K);
                    feasible(q,r) = 1;
                else
                    pok = pk(:,k-1,n);
                    vok = vk(:,k-1,n);
                    aok = ak(:,k-1,n);
                    [pi,vi,ai,feasible(q,r),outbound(q,r),coll(q,r)] = solveSoftDMPCbound(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,A_p_dmpc,A_v_dmpc,Delta,Q,S,E1,E2,order,term); 
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
                fail = fail + 1;
                break;
            end
            l = new_l;
            reached_goal = ReachedGoal(pk,pf,k,error_tol,N);
            k = k+1;
        end
        
        if feasible(q,r) && reached_goal
            at_goal = 1;
        elseif feasible(q,r) && ~reached_goal
            failed_goal(q,r) = 1;
        end

        if feasible(q,r) && ~failed_goal(q,r) 
            % scale the trajectory to meet the limits and plot
            vmax = 2;
            amax = 1;
            ak_mod = [];
            vk_mod = [];
            for i=1:N
                ak_mod(:,i) = amax./sqrt(sum(ak(:,:,i).^2,1));
                vk_mod(:,i) = vmax./sqrt(sum(vk(:,:,i).^2,1));
            end
            r_factor = min([min(min(ak_mod)), min(min(vk_mod))]);
            h_scaled = h/sqrt(r_factor);

            % Time settings and variables
            T = (k-2)*h_scaled; % Trajectory final time
            tk = 0:h_scaled:T;
            Ts = 0.01; % period for interpolation @ 100Hz
            t = 0:Ts:T; % interpolated time vector
            K = T/h_scaled + 1;
            
            % Compute new velocity and acceleration profiles
            for i = 1:N
                for k = 1:size(pk,2)-1
                    ak(:,k,i) = ak(:,k,i)*r_factor;
                    vk(:,k+1,i) = vk(:,k,i) + h_scaled*ak(:,k,i);
                    pk(:,k+1,i) = pk(:,k,i) + h_scaled*vk(:,k,i) + h_scaled^2/2*ak(:,k,i);
                end
            end

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
    end
end
fprintf("Finished! \n")
save('failure_rate_newctr')
%% Post-Processing
close all

% If for any q,r coll and violation =1, then a collision was detected and
% it was not actually infeasible... Make violation = 0
% 
% for q = 1:length(N_vector)
%     for r = 1:trials
%         if (coll(q,r)==1 && feasible(q,r)==0)
%             feasible(q,r)=1; 
%         end
%     end
% end

% Probability of success plots
prob_dmpc = sum(success_dmpc,2)/trials;

figure(1)
grid on;
hold on;
ylim([0,1.05])
plot(N_vector,prob_dmpc,'Linewidth',2);
xlabel('Number of agents');
ylabel('Success Probability');

% Computation time
tmean_dmpc = nanmean(t_dmpc,2);
tstd_dmpc = nanstd(t_dmpc,1,2);
figure(2)
grid on;
hold on;
errorbar(N_vector,tmean_dmpc,tstd_dmpc,'Linewidth',2);
xlabel('Number of agents');
ylabel('Average Computation time [s]');

% Completion time
tmean_traj = nanmean(traj_time,2);
tstd_traj = nanstd(traj_time,1,2);

figure(3)
grid on;
hold on;
errorbar(N_vector,tmean_traj,tstd_traj,'Linewidth',2);
xlabel('Number of agents');
ylabel('Average Time for Transition [s]');

% Failure analysis
violation_num = sum(violation,2);
goal_num = sum(failed_goal,2);
infes_num = sum(~feasible,2);
total_num = sum(violation_num) + sum(goal_num) + sum(infes_num);

StackData = [infes_num/trials*100 violation_num/trials*100 goal_num/trials*100];

figure(4)
h=bar(N_vector,StackData,'stacked','BarWidth',0.3);
myC= summer(size(StackData,2));
for k = 1:size(StackData,2)
    set(h(k),'facecolor',myC(k,:));
end

xticks(N_vector);
set(gca,'fontsize',14)
% set(gca,'LineWidth',2,'TickLength',[0.015 0.015]);
xlabel('Number of agents');
ylabel('Failure Rate %');
legend('Infeasibility','Collisions','Incomplete Trajectory')
% Print various statistics
fprintf("--------------------------- \n k_ctr = k \n");
fprintf("--------------------------- \n")
fprintf("Probability of success across all tests --> %.2f%% (%d out of %d tests) \n",...
    sum(prob_dmpc)/length(prob_dmpc)*100, sum(sum(success_dmpc)), length(N_vector)*trials)
fprintf("Percentage failure due to infeasibility --> %.2f%% (%d out of %d failures) \n",...
    sum(infes_num)/total_num*100, sum(infes_num), total_num)

fprintf("Percentage failure due to collisions --> %.2f%% (%d out of %d failures) \n",...
    sum(violation_num)/total_num*100, sum(violation_num), total_num)

fprintf("Percentage failure due to not reaching goal --> %.2f%% (%d out of %d failures) \n",...
    sum(goal_num)/total_num*100, sum(goal_num), total_num)
