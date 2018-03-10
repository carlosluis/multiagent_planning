function [p,v,a] = singleSCP(po, pf, h, K, B, l)

p = initSolution(po,pf,h,K);
epsilon = 2; %to be tuned
criteria = 0;

while criteria > epsilon
    
    % Setup the QP
    [Ain, bin] = IneqConstr(p, l, h);
    [Aeq, beq] = EqConstr(pf,K,h);
    H = eye(3*K);
    
    
    
end

end