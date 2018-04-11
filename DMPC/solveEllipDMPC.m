function [p,v,a,success] = solveEllipDMPC(po,pf,vo,ao,n,h,l,K,rmin,pmin,pmax,alim,A,A_initp,Delta,tol,Q1,S1)

success = 1;
k_hor = size(l,2);
val = tol + 2;
ub = alim*ones(3*K,1);
lb = -ub; 
i = 1;
c = 1.5;
E = diag([1,1,c]);
addConstr = [];
prev_p = l(:,:,n);
Aeq = [];
beq = [];
Ain_total = [];
bin_total = [];
options = optimset('Display', 'off');

while (i <= k_hor && val > tol)
    newConstrCount = 0; 
    Ain_total = [];
    bin_total = [];
    for k = 1: k_hor
        violation = CheckCollEllipDMPC(prev_p(:,k),l,n,k,E);
        if (ismember(k,addConstr))
            [Ainr, binr] = CollConstrEllipDMPC(prev_p(:,k),po,vo,n,k,l,A,E,A_initp);
            Ain_total = [Ain_total; Ainr];
            bin_total = [bin_total; binr];
            
        elseif (newConstrCount==0 && violation)
            [Ainr, binr] = CollConstrEllipDMPC(prev_p(:,k),po,vo,n,k,l,A,E,A_initp);
            Ain_total = [Ain_total; Ainr];
            bin_total = [bin_total; binr];  
            addConstr = [addConstr k];
            newConstrCount = newConstrCount + 1;
        end       
    end
    
    % Setup the QP
    if(isempty(Ain_total)) % Case of no collisions
        Q = 1000*[zeros(3*(K-1),3*K);
                zeros(3,3*(K-1)) eye(3)];
        R = 1*eye(3*K);
        S = 10*eye(3*K);
    else
        Q = Q1*[zeros(3*(K-1),3*K);
                zeros(3,3*(K-1)) eye(3)];
        R = 1*eye(3*K);
        S = S1*eye(3*K);
    end
    
    Ain_total = [Ain_total; A; -A];
    bin_total = [bin_total; repmat((pmax)',K,1) - A_initp*([po';vo']); repmat(-(pmin)',K,1) + A_initp*([po';vo'])];
    H = 2*(A'*Q*A+ Delta'*S*Delta + R);
    ao_1 = [ao zeros(1,3*(K-1))];
    f = -2*(repmat((pf)',K,1)'*Q*A - (A_initp*([po';vo']))'*Q*A + ao_1*S*Delta) ;
    
    %Solve and propagate states
    [a,fval,exitflag] = quadprog(H,f',Ain_total,bin_total,Aeq,beq,lb,ub,[],options);
    if (isempty(a) || exitflag == 0)
        p = [];
        v = [];
        success = 0;
        return
    end
    success = exitflag;
    [p,v] = propStatedmpc(po,vo,a,h);
    p = vec2mat(p,3)';
    v = vec2mat(v,3)';
    a = vec2mat(a,3)';
    val = maxDeviation(p, prev_p);  
    prev_p = p;
    i = i + 1;     
end
% fprintf("Number of SCP iterations = %i\n",i-1)
end


