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
success = 1;
N_vector = 2:2:30; % number of vehicles
trials = 50;
fail = 0;
% Workspace boundaries
pmin = [-2.5,-2.5,0.2];
pmax = [2.5,2.5,2.2];

% Minimum distance between vehicles in m
rmin = 0.75;

% Maximum acceleration in m/s^2
alim = 0.7;

% Some Precomputations dec-iSCP
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
A_p = [];
A_v = [];

% Build matrix to convert acceleration to position
for k = 1:(K-1)
    add_b = [zeros(size(b,1),size(b,2)*(k-1)) b zeros(size(b,1),size(b,2)*(K-k))];
    new_row = A*prev_row + add_b;   
    A_p = [A_p; new_row(1:3,:)];
    A_v = [A_v; new_row(4:6,:)];
    prev_row = new_row; 
end

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

        % Empty list of obstacles
        l = [];
        
        % DEC-ISCP
        t_start = tic; 
        for i = 1:N 
            poi = po(:,:,i);
            pfi = pf(:,:,i);
            [pi, vi, ai,success] = singleiSCP(poi,pfi,h,K,pmin,pmax,rmin,alim,l,A_p,A_v);
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
        tries(q,r) = 1;
        failed_goal(q,r) = 0;
        t_start = tic;
        while tries(q,r) <= 10 && ~at_goal
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
                    tries(q,r) = tries(q,r) + 1;
                    Q = Q+50;
                    break;
                end
                l = new_l;
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
            for i = 1:N
                p(:,:,i) = spline(tk,pk(:,:,i),t);
                v(:,:,i) = spline(tk,vk(:,:,i),t);
                a(:,:,i) = spline(tk,ak(:,:,i),t); 
            end
            t_dmpc(q,r) = toc(t_start);
            totdist_dmpc(q,r) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
        else
            t_dmpc(q,r) = nan;
            totdist_dmpc(q,r) = nan;
            save(['Fail_' num2str(fail)]);
            fail = fail + 1;
        end
        success_dmpc(q,r) = success && at_goal;
    end
end
fprintf("Finished!")
%% Post-Processing

% Probability of success plots
prob_dec = sum(success_dec,2)/trials;
prob_dmpc = sum(success_dmpc,2)/trials;
figure(1)
plot(N_vector,prob_dec','Linewidth',2);
grid on;
hold on;
ylim([0,1.05])
xlim([0,27])
plot(N_vector,prob_dmpc,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Success Probability');
legend('dec-iSCP','DMPC');

% Computation time
tmean_dec = nanmean(t_dec,2);
tstd_dec = nanstd(t_dec,1,2);
tmean_dmpc = nanmean(t_dmpc,2);
tstd_dmpc = nanstd(t_dmpc,1,2);
figure(2)
errorbar(N_vector,tmean_dec,tstd_dec,'Linewidth',2);
grid on;
hold on;
xlim([0,27])
errorbar(N_vector,tmean_dmpc,tstd_dmpc,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Average Computation Time [s]');
legend('dec-iSCP','DMPC');

% Percentage increase/decrease on travelled dist of dmpc wrt dec
% Positive number means that dmpc path was longer
diff_dist = (totdist_dmpc-totdist_dec)./totdist_dmpc;
avg_diff = nanmean(diff_dist,2);
std_diff = nanstd(diff_dist,1,2);
figure(3)
errorbar(N_vector,100*avg_diff,100*std_diff,'Linewidth',2);
grid on;
xlabel('Number of Vehicles');
ylabel('Average % increase/decrease');
title('Percentual increase/decrease on total travelled distance of DMPC wrt dec-iSCP');

