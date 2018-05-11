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
trials = 30;
% Workspace boundaries
pmin = [-2.5,-2.5,0.2];
pmax = [2.5,2.5,2.2];

% Minimum distance between vehicles in m
rmin = 0.75;

% Maximum acceleration in m/s^2
alim = 0.7;

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

Q = [1000,1000,100,1000,10000];
S = [10,100,100,1000,10];
tol = 2;

for q = 1:length(N_vector)
    N = N_vector(q);
    for r = 1:trials
        fprintf("Doing trial #%i with %i vehicles\n",r,N)
        % Initial positions
        [po,pf] = randomTest(N,pmin,pmax,rmin);
        for m = 1:length(Q)
            %DMPC just one iteration i.e really high tolerance
            % Empty list of obstacles
            l = [];
            t_start = tic;
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
                        [pi,vi,ai,success] = solveDMPC(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta,tol,Q(m),S(m)); 
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

            if success
                for i = 1:N
                    p(:,:,i) = spline(tk,pk(:,:,i),t);
                    v(:,:,i) = spline(tk,vk(:,:,i),t);
                    a(:,:,i) = spline(tk,ak(:,:,i),t); 
                end
                t_dmpc(q,r,m) = toc(t_start);
                totdist_dmpc(q,r,m) = sum(sum(sqrt(diff(p(1,:,:)).^2+diff(p(2,:,:)).^2+diff(p(3,:,:)).^2)));
            else
                t_dmpc(q,r,m) = nan;
                totdist_dmpc(q,r,m) = nan;
            end
            success_dmpc(q,r,m) = success;
        end     
    end
end
fprintf("Finished!\n")
%% Post-Processing

% Probability of success plots
figure(1)
for i = 1:length(Q)
    prob_dmpc(:,i) = sum(success_dmpc(:,:,i),2)/trials;
    h_plot(i) = plot(N_vector,prob_dmpc(:,i),'Linewidth',2);
    h_label{i} = ['Q = ' num2str(Q(i)) ', S = ' num2str(S(i))];
    grid on;
    hold on;
end
ylim([0,1.05])
xlabel('Number of Vehicles');
ylabel('Success Probability');
legend(h_plot,h_label);

% Computation time
figure(2)
for i = 1:length(Q)
    tmean_dmpc = nanmean(t_dmpc(:,:,i),2);
    h_plot(i) = plot(N_vector,tmean_dmpc,'Linewidth',2);
    grid on;
    hold on;
end
xlabel('Number of Vehicles');
ylabel('Average Computation time [s]');
legend(h_plot,h_label);