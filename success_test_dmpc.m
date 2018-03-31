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
fail = 0;

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
        
        %DMPC
        % Empty list of obstacles
        l = [];
        success = 0;
        Q = 100;
        S = 100; 
        tries(q,r) = 1;
        t_start = tic;
        while tries(q,r) <= 10 && ~success
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
                    Q = Q+100;
                    if (tries > 5)
                        S = S - 20;
                    end
                    break;
                end
                l = new_l;
            end   
        end

        if success
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
        success_dmpc(q,r) = success;
    end
end
fprintf("Finished!")
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
legend('DMPC');

% Computation time
tmean_dmpc = nanmean(t_dmpc,2);
tstd_dmpc = nanstd(t_dmpc,1,2);
figure(2)
grid on;
hold on;
errorbar(N_vector,tmean_dmpc,tstd_dmpc,'Linewidth',2);
xlabel('Number of Vehicles');
ylabel('Average Computation time [s]');
legend('DMPC');
