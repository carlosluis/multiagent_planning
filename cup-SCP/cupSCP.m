clc
clear all
close all

% Time settings and variables
T = 15; % Trajectory final time
h = 0.2; % time step duration
tk = 0:h:T;
K = T/h + 1; % number of time steps
Ts = 0.01; % period for interpolation @ 100Hz
t = 0:Ts:T; % interpolated time vector
success = 1;
N = 2; % number of vehicles

% Workspace boundaries
pmin = [-2.5,-2.5,0.2];
pmax = [2.5,2.5,2.2];

% Minimum distance between vehicles in m
rmin = 0.75;

% Maximum acceleration in m/s^2
alim = 0.5;

% Initial positions
[po,pf] = randomTest(N,pmin,pmax,rmin);

[p,v,a] = solveCupSCP(po,pf,h,K,N,pmin,pmax,rmin,alim);
