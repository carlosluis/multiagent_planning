function [p,v,a] = singleSCP(po, pf, h, K, pmin,pmax, l)

prev_p = initSolution(po,pf,h,K);
epsilon = 2; %to be tuned
tol = 0;
H = eye(3*K);
a_lim = 1; %Maximum acc of 1 m/s^2
ub = a_lim*ones(3*K,1);
lb = -ub; 
i = 1;

while i <= 5
    
    % Setup the QP
    [Ain, bin] = IneqConstr(prev_p, l, h, pmin, pmax);
    [Aeq, beq] = EqConstr(pf-po,K,h);
    
    %Solve and propagate states
    a = quadprog(H,[],Ain,bin,Aeq,beq,lb,ub);   
    [p,v] = propState(po,a,h);
    p = vec2mat(p,3)';
    v = vec2mat(v,3)';
    a = vec2mat(a,3)';
    tol = maxDeviation(p, prev_p);
    
    prev_p = p;
    i = i + 1;
    
end

end