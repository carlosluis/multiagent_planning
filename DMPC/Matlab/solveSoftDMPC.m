function [p,v,a,success,outbound] = solveSoftDMPC(po,pf,vo,ao,n,h,l,K,rmin,pmin,pmax,alim,A,A_initp,Delta,Q1,S1,E1,E2,order)

k_hor = size(l,2);
ub = alim*ones(3*K,1);
lb = -ub; 
prev_p = l(:,:,n);
Aeq = [];
beq = [];
options =  optimset('Display', 'off');
N = size(l,3);
Ain_coll = []; 
bin_coll = []; 

for k = 1: k_hor
    violation = CheckCollEllipDMPC(prev_p(:,k),l,n,k,E1,rmin,order);
    if (violation)
        [Ainr, binr] = CollConstrEllipDMPC(prev_p(:,k),po,vo,n,k,l,rmin,A,A_initp,E1,E2,order);
        Ain_coll = [Ainr eye(N-1,N-1);
                     zeros(N-1,3*K) eye(N-1)];
        bin_coll = [binr; zeros(N-1,1)];
        break;
    end       
end

spd = 1;

% Setup the QP
if(isempty(Ain_coll) && norm(po-pf) >= 1) % Case of no collisions far from sp
    Q = 1000*[zeros(3*(K-spd),3*K);
            zeros(3*spd,3*(K-spd)) eye(3*spd)];
    R = 1*eye(3*K);
    S = 10*eye(3*K);
elseif (isempty(Ain_coll) && norm(po-pf) < 1) % no collisions close to sp
    Q = 10000*[zeros(3*(K-spd),3*K);
            zeros(3*spd,3*(K-spd)) eye(3*spd)];
    R = 1*eye(3*K);
    S = 10*eye(3*K); 
else     % collisions
    Q = Q1*[zeros(3*(K-spd),3*K);
            zeros(3*spd,3*(K-spd)) eye(3*spd)];
    R = 1*eye(3*K);
    S = S1*eye(3*K);
end

if (violation) % In case of collisions, we relax the constraint with slack variable
    % Add dimensions for slack variable
    Q = [Q zeros(3*K,N-1);
         zeros(N-1,3*K) zeros(N-1,N-1)];
    R = [R zeros(3*K,N-1);
         zeros(N-1,3*K) zeros(N-1,N-1)];
    S = [S zeros(3*K,N-1);
         zeros(N-1,3*K) zeros(N-1,N-1)];
    A = [A zeros(3*K,N-1);
         zeros(N-1,3*K) zeros(N-1,N-1)];
    Delta = [Delta zeros(3*K,N-1);
         zeros(N-1,3*K) zeros(N-1,N-1)];
    bin_total = [bin_coll; repmat((pmax)',K,1) - A_initp*([po';vo']); zeros(N-1,1); repmat(-(pmin)',K,1) + A_initp*([po';vo']); zeros(N-1,1)];
    ao_1 = [ao zeros(1,3*(K-1)+N-1)];
    A_initp = [A_initp; zeros(N-1,6)];
    
    % Linear penalty on collision constraint relaxation
    f_eps = -5*10^4*[zeros(3*K,1); ones(N-1,1)]';
    
    % Quadratic penalty on collision constraint relaxation
    EPS = 1*10^3*[zeros(3*K,3*K) zeros(3*K,N-1);
           zeros(N-1,3*K) eye(N-1,N-1)];
       
    f = -2*([repmat((pf)',K,1); zeros(N-1,1)]'*Q*A - (A_initp*([po';vo']))'*Q*A + ao_1*S*Delta) + f_eps ;
    
else % case of no collisions, we don't even add collision constraints
    bin_total = [bin_coll; repmat((pmax)',K,1) - A_initp*([po';vo']); repmat(-(pmin)',K,1) + A_initp*([po';vo'])];
    ao_1 = [ao zeros(1,3*(K-1))];
    f = -2*(repmat((pf)',K,1)'*Q*A - (A_initp*([po';vo']))'*Q*A + ao_1*S*Delta) ;
    EPS = zeros(3*K,3*K); 
end

Ain_total = [Ain_coll; A; -A];
H = 2*(A'*Q*A+ Delta'*S*Delta + R + EPS);
outbound = 0;
%Solve and propagate states
[x,fval,exitflag] = quadprog(H,f',Ain_total,bin_total,Aeq,beq,lb,ub,[],options);
if (exitflag == -2 || exitflag == 0 || exitflag == -8) 
    p = [];
    v = [];
    a = [];
    success = 0;
    if (~violation)
        outbound = 1;
    end
    return
end
success = exitflag;
a = x(1:3*K);
if violation % extract the value of the slack variable (not used atm)
    epsilon = x(3*K+1:end);
end
[p,v] = propStatedmpc(po,vo,a,h);
p = vec2mat(p,3)';
v = vec2mat(v,3)';
a = vec2mat(a,3)';     
end