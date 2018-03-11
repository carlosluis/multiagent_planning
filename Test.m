%Testing functions of dec-SCP
T = 2;
h = 0.2;
K = T/h + 1;
po = [0,0,0];
pf = [0,3,4];

pmin = [0,0,0];
pmax = [0,3,4];

p = initSolution(po,pf,h,K);

% plot3(p(1,:), p(2,:), p(3,:) )

p_obs1 = p+1;
p_obs2 = p+2;

p_obs = cat(3, p_obs1, p_obs2, p_obs1, p_obs2, p_obs1, p_obs2);

[Ain, bin] = IneqConstr(p, p_obs, h, pmin, pmax);

[Aeq, beq] = EqConstr(pf,K,h);

H = eye(3*K);

%% Test Single SCP

T = 5;
h = 0.05;
t = 0:h:5;
K = T/h + 1;
po = [0,0,1.5];
pf = [1,2,1];

pmin = [-4,-4,0];
pmax = [4,4,3.5];

l = [];

[p,v,a] = singleSCP(po, pf-po, h, K, pmin,pmax, l);

% Interpolate to obtain more points
Ts = 0.01; % Produce trajectory @ 100Hz
t2 = 0:Ts:5;

p = spline(t,p,t2);
v = spline(t,v,t2);
a = spline(t,a,t2);

figure(1)
plot3(p(1,:), p(2,:), p(3,:));
grid on;

figure(2)
subplot(3,1,1)
plot(t2,p(1,:));
ylabel('x [m]')
xlabel ('t [s]')
grid on;

subplot(3,1,2)
plot(t2,p(2,:));
ylabel('y [m]')
xlabel ('t [s]')
grid on;

subplot(3,1,3)
plot(t2,p(3,:));
ylabel('z [m]')
xlabel ('t [s]')
grid on;

figure(3)
subplot(3,1,1)
plot(t2,v(1,:));
ylabel('vx [m/s]')
xlabel ('t [s]')
grid on;

subplot(3,1,2)
plot(t2,v(2,:));
ylabel('vy [m/s]')
xlabel ('t [s]')
grid on;

subplot(3,1,3)
plot(t2,v(3,:));
ylabel('vz [m/s]')
xlabel ('t [s]')
grid on;

figure(4)
subplot(3,1,1)
plot(t2,a(1,:));
ylabel('ax [m/s]')
xlabel ('t [s]')
grid on;

subplot(3,1,2)
plot(t2,a(2,:));
ylabel('ay [m/s]')
xlabel ('t [s]')
grid on;

subplot(3,1,3)
plot(t2,a(3,:));
ylabel('az [m/s]')
xlabel ('t [s]')
grid on;








