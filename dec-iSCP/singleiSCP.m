function [p,v,a,success] = singleiSCP(po,pf,h,K,pmin,pmax,rmin,alim,l,A_p,A_v,E1,E2,order)

prev_p = initSolution(po,pf,h,K);
ub = alim*ones(3*K,1);
lb = -ub; 
i = 1;
success = 1;
H = eye(3*K);
A = getPosMat(h,K);
Aeq = getPosVelMat(h,K);
options = optimoptions('quadprog','Display','off');
check = true;
addConstr = [];

while (i <= K && check)
    newConstrCount = 0; 
    Ain_total = [];
    bin_total = [];
    for k = 1:K
        violation = CheckforColl(prev_p(:,k),l,k,E1,rmin,order);
        if (ismember(k,addConstr))
            [Ainr, binr] = CollConstr(prev_p(:,k),po,k,l,A,rmin,E1,E2,order);
            Ain_total = [Ain_total; Ainr];
            bin_total = [bin_total; binr];
        
        elseif (newConstrCount==0 && violation)
            [Ainr, binr] = CollConstr(prev_p(:,k),po,k,l,A,rmin,E1,E2,order);
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
    [a,fval,exitflag] = quadprog(H,[],Ain_total,bin_total,Aeq,beq,lb,ub,[],options);
    if (isempty(a))
        p = [];
        v = [];
        success = 0;
        return
    end
    [p,v] = propState(po,a,A_p,A_v,K);
    p = vec2mat(p,3)';
    v = vec2mat(v,3)';
    a = vec2mat(a,3)';
    fail = CheckforAllColl(p,l,K,E1,rmin,addConstr,order);
    if ~fail
        check = false;
    end
    prev_p = p;
    i = i + 1;   
end
if (check)
    success = 0;
end
% fprintf("Number of SCP iterations = %i\n",i-1)
end