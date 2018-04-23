function [p,v,a] = singleSCP(po,pf,h,K,pmin,pmax,rmin,alim,l)
po = squeeze(po);
prev_p = initSolution(po,pf,h,K);
epsilon = 2; %to be tuned
tol = 2;
ub = alim*ones(3*K,1);
lb = -ub; 
i = 1;

H = eye(3*K);
A = getPosMat(h,K);
Aeq = getPosVelMat(h,K);

while (i <= K && tol > 0.01)
    
    % Setup the QP
    [Ain, bin] = AllCollConstr(prev_p,K,rmin,l,A);
    Ain_total = [Ain; A; -A];
    bin_total = [bin; repmat((pmax-po)',K,1); repmat(-(pmin-po)',K,1)];
    beq = [(pf-po)' ; zeros(3,1); zeros(3,1); zeros(3,1)];
    
    %Solve and propagate states
    a = quadprog(H,[],Ain_total,bin_total,Aeq,beq,lb,ub);   
    [p,v] = propState(po,a,h);
    p = vec2mat(p,3)';
    v = vec2mat(v,3)';
    a = vec2mat(a,3)';
    tol = maxDeviation(p, prev_p);
    
    prev_p = p;
    i = i + 1;    
end

end