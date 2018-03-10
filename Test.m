%Testing functions of dec-SCP
T = 2;
h = 0.2;
K = T/h + 1;
po = [0,0,0];
pf = [0,3,4];

pmin = [0,0,0];
pmax = [0,3,4];

p = initSolution(po,pf,h,K);

plot3(p(1,:), p(2,:), p(3,:) )

p_obs1 = p+1;
p_obs2 = p+2;

p_obs = cat(3, p_obs1, p_obs2, p_obs1, p_obs2, p_obs1, p_obs2);

[Ain, bin] = IneqConstr(p, p_obs, h, pmin, pmax);

[Aeq, beq] = EqConstr(pf,K,h);

H = eye(3*K);




