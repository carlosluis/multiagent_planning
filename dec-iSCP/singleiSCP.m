function [p,v,a] = singleiSCP(po,pf,h,K,pmin,pmax,rmin,alim,l)

prev_p = initSolution(po,pf,h,K);
tol = 2;
ub = alim*ones(3*K,1);
lb = -ub; 
i = 1;

H = eye(3*K);
A = getPosMat(h,K);
Aeq = getPosVelMat(h,K);

addConstr = [];

while (i <= K && tol > 0.001)
    newConstrCount = 0; 
    Ain_total = [];
    bin_total = [];
    for k = 1:K
        violation = CheckforColl(prev_p(:,k),l,k,rmin);
        if (ismember(k,addConstr))
            [Ainr, binr] = CollConstr(prev_p(:,k),po,k,l,A,rmin);
            Ain_total = [Ain_total; Ainr];
            bin_total = [bin_total; binr];
        
        elseif (newConstrCount==0 && violation)
            [Ainr, binr] = CollConstr(prev_p(:,k-1),po,k,l,A,rmin);
            Ain_total = [Ain_total; Ainr];
            bin_total = [bin_total; binr];  
            addConstr = [addConstr k];
            newConstrCount = newConstrCount + 1;
        end           
    end
    
    % Setup the QP
    Ain_total = [Ain_total; A; -A];
    bin_total = [bin_total; repmat((pmax-po)',K,1); repmat(-(pmin-po)',K,1)];
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
fprintf("Number of SCP iterations = %i\n",i-1)
end