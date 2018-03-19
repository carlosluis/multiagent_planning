function [p,v,a] = solveDMPC(po,pf,vo,ao,n,h,k_step,l,K,rmin,pmin,pmax,alim)

k_hor = size(l,2);
A = getPosMat(h,K);
tol = 2;
ub = alim*ones(3*K,1);
lb = -ub; 
i = 1;
addConstr = [];
prev_p = l(:,:,n);
Aeq = [eye(3) zeros(3,3*(K-1))];
beq = ao';

while (i <= k_hor && tol > 0.01)
    newConstrCount = 0; 
    Ain_total = [];
    bin_total = [];
    for k = 1: k_hor
        violation = CheckCollDMPC(prev_p,l,n,k,rmin);
        if (ismember(k,addConstr))
            [Ainr, binr] = CollConstrDMPC(prev_p(:,k),po,n,k,l,A,rmin);
            Ain_total = [Ain_total; Ainr];
            bin_total = [bin_total; binr];
            
        elseif (newConstrCount==0 && violation)
            [Ainr, binr] = CollConstrDMPC(prev_p(:,k),po,n,k,l,A,rmin);
            Ain_total = [Ain_total; Ainr];
            bin_total = [bin_total; binr];  
            addConstr = [addConstr k];
            newConstrCount = newConstrCount + 1;
        end       
    end
    
    % Setup the QP
    Ain_total = [Ain_total; A; -A];
    bin_total = [bin_total; repmat(pmax',K,1); repmat(-pmin',K,1)];
    Q = eye(3*K);
    R = eye(3*K);
    H = 2*(A'*Q*A);
    f = -2*repmat(pf',K,1)'*Q*A;
    
    %Solve and propagate states
    a = quadprog(H,f',Ain_total,bin_total,Aeq,beq,lb,ub);   
    [p,v] = propStatedmpc(po,vo,a,h);
    p = vec2mat(p,3)';
    v = vec2mat(v,3)';
    a = vec2mat(a,3)';
    tol = maxDeviation(p, prev_p);  
    prev_p = p;
    i = i + 1;     
end


