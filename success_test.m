clc
clear all
close all

% Time settings and variables
T = 12; % Trajectory final time
h = 0.2; % time step duration
tk = 0:h:T;
K = T/h + 1; % number of time steps
Ts = 0.01; % period for interpolation @ 100Hz
t = 0:Ts:T; % interpolated time vector
k_hor = 15;
success = 1;
N_vector = 2; % number of vehicles
trials = 5;
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

for q = 1:length(N_vector)
    N = N_vector(q);
    for r = 1:trials
        % Initial positions
        [po,pf] = randomTest(N,pmin,pmax,rmin);

        % Empty list of obstacles
        l = [];
        
        % DEC-ISCP
        t_start = tic; 
        for i = 1:N 
            poi = po(:,:,i);
            pfi = pf(:,:,i);
            [pi, vi, ai,success] = singleiSCP(poi,pfi,h,K,pmin,pmax,rmin,alim,l);
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
            t_dec(q,r) = 0;
            totdist_dec(q,r) = 0;
        end
        success_dec(q,r) = success;
        
        %DMPC
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
                    [pi,vi,ai,success] = solveDMPC(pok',pf(:,:,n),vok',aok',n,h,l,k_hor,rmin,pmin,pmax,alim,A,A_initp,Delta); 
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
            pred(:,:,:,k) = l;
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
            t_dmpc(q,r) = 0;
            totdist_dmpc(q,r) = 0;
        end
        success_dmpc(q,r) = success;
    end
end