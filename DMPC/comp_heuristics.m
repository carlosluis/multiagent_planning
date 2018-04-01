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
k_hor = 15;
N_vector = 2:2:30; % number of vehicles
trials = 50;
fail2 = 0;
fail3 = 0;
% Workspace boundaries
pmin = [-2.5,-2.5,0.2];
pmax = [2.5,2.5,2.2];

% Minimum distance between vehicles in m
rmin = 0.75;

% Maximum acceleration in m/s^2
alim = 0.7;

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
        [po,pf] = randomTest(N,pmin,pmax,rmin);
        
        %DMPC No heuristic
        % Empty list of obstacles
        l = [];
        success = 0;
        epsilon = 0;
        Q = 1000;
        S = 10; 
        t_start = tic;
        for k = 1:K
            for n = 1:N
                if k==1
                    poi = po(:,:,n);
                    pfi = pf(:,:,n);
                    [pi,vi,ai] = initDMPC(poi,pfi,h,k_hor,K,epsilon);
                    success = 1;
                else
                    pok = pk(:,k-1,n);
                    vok = vk(:,k-1,n);
                    aok = ak(:,k-1,n);
                    [pi,vi,ai,success] = solveDMPC(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta,tol,Q,S); 
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
                break;
            end
            l = new_l;  
        end

        if success && ReachedGoal(pk,pf,K,error_tol)
            for i = 1:N
                p(:,:,i) = spline(tk,pk(:,:,i),t);
                v(:,:,i) = spline(tk,vk(:,:,i),t);
                a(:,:,i) = spline(tk,ak(:,:,i),t); 
            end
            t_dmpc1(q,r) = toc(t_start);
            totdist_dmpc1(q,r) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
        else
            t_dmpc1(q,r) = nan;
            totdist_dmpc1(q,r) = nan;
        end
        success_dmpc1(q,r) = success && ReachedGoal(pk,pf,K,error_tol);
        
        %DMPC No epsilon heuristic
        % Empty list of obstacles
        tol = 2;
        success = 0; %check if QP was feasible
        at_goal = 0; %At the end of solving, makes sure every agent arrives at the goal
        error_tol = 0.05; % 5cm destination tolerance
        epsilon = 0; % heuristic variable to initialize DMPC more conservative

        % Penalty matrices when there're predicted collisions
        Q = 10;
        S = 100; 
        tries2(q,r) = 1;
        failed_goal2(q,r) = 0;
        t_start = tic;
        while tries2(q,r) <= 10 && ~at_goal
            for k = 1:K
                for n = 1:N
                    if k==1
                        poi = po(:,:,n);
                        pfi = pf(:,:,n);
                        [pi,vi,ai] = initDMPC(poi,pfi,h,k_hor,K,epsilon);
                        success = 1;
                    else
                        pok = pk(:,k-1,n);
                        vok = vk(:,k-1,n);
                        aok = ak(:,k-1,n);
                        [pi,vi,ai,success] = solveDMPC(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta,tol,Q,S); 
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
                    tries2(q,r) = tries2(q,r) + 1;
                    Q = Q+50;
%                     epsilon = epsilon + 5;
                    break;
                end
                l = new_l;
            end   
            pass = ReachedGoal(pk,pf,K,error_tol);
            if success && pass
                at_goal = 1;
            elseif success && ~pass
                failed_goal2(q,r) = failed_goal2(q,r) + 1;
                tries2(q,r) = tries2(q,r) + 1;
                Q = Q+100;
            end
        end

        if success && at_goal
            for i = 1:N
                p(:,:,i) = spline(tk,pk(:,:,i),t);
                v(:,:,i) = spline(tk,vk(:,:,i),t);
                a(:,:,i) = spline(tk,ak(:,:,i),t); 
            end
            t_dmpc2(q,r) = toc(t_start);
            totdist_dmpc2(q,r) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
        else
            t_dmpc2(q,r) = nan;
            totdist_dmpc2(q,r) = nan;
            save(['Fail2_' num2str(fail2)]);
            fail2 = fail2 + 1;
        end
        success_dmpc2(q,r) = success && at_goal;
        
        %DMPC All heuristics
        % Empty list of obstacles
        tol = 2;
        success = 0; %check if QP was feasible
        at_goal = 0; %At the end of solving, makes sure every agent arrives at the goal
        error_tol = 0.05; % 5cm destination tolerance
        epsilon = 0; % heuristic variable to initialize DMPC more conservative

        % Penalty matrices when there're predicted collisions
        Q = 10;
        S = 100; 
        tries3(q,r) = 1;
        failed_goal3(q,r) = 0;
        t_start = tic;
        while tries3(q,r) <= 10 && ~at_goal
            for k = 1:K
                for n = 1:N
                    if k==1
                        poi = po(:,:,n);
                        pfi = pf(:,:,n);
                        [pi,vi,ai] = initDMPC(poi,pfi,h,k_hor,K,epsilon);
                        success = 1;
                    else
                        pok = pk(:,k-1,n);
                        vok = vk(:,k-1,n);
                        aok = ak(:,k-1,n);
                        [pi,vi,ai,success] = solveDMPC(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta,tol,Q,S); 
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
                    tries3(q,r) = tries3(q,r) + 1;
                    Q = Q+50;
                    epsilon = epsilon + 5;
                    break;
                end
                l = new_l;
            end   
            pass = ReachedGoal(pk,pf,K,error_tol);
            if success && pass
                at_goal = 1;
            elseif success && ~pass
                failed_goal3(q,r) = failed_goal3(q,r) + 1;
                tries3(q,r) = tries3(q,r) + 1;
                Q = Q+100;
            end
        end

        if success && at_goal
            for i = 1:N
                p(:,:,i) = spline(tk,pk(:,:,i),t);
                v(:,:,i) = spline(tk,vk(:,:,i),t);
                a(:,:,i) = spline(tk,ak(:,:,i),t); 
            end
            t_dmpc3(q,r) = toc(t_start);
            totdist_dmpc3(q,r) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
        else
            t_dmpc3(q,r) = nan;
            totdist_dmpc3(q,r) = nan;
            save(['Fail3_' num2str(fail3)]);
            fail3 = fail3 + 1;
        end
        success_dmpc3(q,r) = success && at_goal;
    end
end
fprintf("Finished!")
%% Post-Processing

% Probability of success plots
prob_dmpc1 = sum(success_dmpc1,2)/trials;
prob_dmpc2 = sum(success_dmpc2,2)/trials;
prob_dmpc3 = sum(success_dmpc3,2)/trials;
figure(1)
grid on;
hold on;
ylim([0,1.05])
plot(N_vector,prob_dmpc1,'Linewidth',2);
plot(N_vector,prob_dmpc2,'Linewidth',2);
plot(N_vector,prob_dmpc3,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Success Probability');
legend('No Heuristics','No Epsilon','All Heuristics');

% Computation time
tmean_dmpc1 = nanmean(t_dmpc1,2);
tstd_dmpc1 = nanstd(t_dmpc1,1,2);
tmean_dmpc2 = nanmean(t_dmpc2,2);
tstd_dmpc2 = nanstd(t_dmpc2,1,2);
tmean_dmpc3 = nanmean(t_dmpc3,2);
tstd_dmpc3 = nanstd(t_dmpc3,1,2);
figure(2)
grid on;
hold on;
errorbar(N_vector,tmean_dmpc1,tstd_dmpc1,'Linewidth',2);
errorbar(N_vector,tmean_dmpc2,tstd_dmpc2,'Linewidth',2);
errorbar(N_vector,tmean_dmpc3,tstd_dmpc3,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Average Computation time [s]');
legend('No Heuristics','No Epsilon','All Heuristics');

% Average number of tries
tries_mean2 = mean(tries2,2);
tries_std2 = std(tries2,1,2);
tries_mean3 = mean(tries3,2);
tries_std3 = std(tries3,1,2);
figure(3)
grid on;
hold on;
errorbar(N_vector,tries_mean2,tries_std2,'LineWidth',2);
errorbar(N_vector,tries_mean3,tries_std3,'LineWidth',2);
xlabel('Number of Vehicles');
ylabel('Average number of DMPC iterations');
legend('No Epsilon','All Heuristics');
