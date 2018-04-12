function [p,v,a,success] = solveEllipDMPC(po,pf,vo,ao,n,h,l,K,rmin,pmin,pmax,alim,A,A_initp,Delta,Q1,S1,E1,E2,order)

k_hor = size(l,2);
ub = alim*ones(3*K,1);
lb = -ub; 
addConstr = [];
prev_p = l(:,:,n);
Aeq = [];
beq = [];
Ain_total = [];
bin_total = [];
options = optimset('Display', 'off');
 
for k = 1: k_hor
    violation = CheckCollEllipDMPC(prev_p(:,k),l,n,k,E1,rmin,order);
    if (violation)
        [Ainr, binr] = CollConstrEllipDMPC(prev_p(:,k),po,vo,n,k,l,rmin,A,A_initp,E1,E2,order);
        Ain_total = [Ain_total; Ainr];
        bin_total = [bin_total; binr];
        break;
    end       
end

% Setup the QP
if(isempty(Ain_total) && norm(po-pf) > 1) % Case of no collisions far from sp
    Q = 100*[zeros(3*(K-1),3*K);
            zeros(3,3*(K-1)) eye(3)];
    R = 1*eye(3*K);
    S = 100*eye(3*K);
elseif (isempty(Ain_total) && norm(po-pf) < 1) % no collisions close to sp
    Q = 1000*[zeros(3*(K-1),3*K);
            zeros(3,3*(K-1)) eye(3)];
    R = 1*eye(3*K);
    S = 1*eye(3*K); 
else     % collisions
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
end


