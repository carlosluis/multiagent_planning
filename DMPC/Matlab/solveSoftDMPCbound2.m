function [p,v,a,success,outbound,coll] = solveSoftDMPCbound2(po,pf,vo,ao,n,h,l,K,rmin,pmin,pmax,alim,A,A_initp,A_p,A_v,Delta,Q1,S1,E1,E2,order,term)

k_hor = size(l,2);
ub = alim*ones(3*K,1);
lb = -ub; 
prev_p = l(:,:,n);
% clip prev_p to within the boundaries
% prev_p = bsxfun(@min,prev_p,pmax');
% prev_p = bsxfun(@max,prev_p,pmin');
constr_tol = 1e-3;
Aeq = [];
beq = [];
options = optimoptions('quadprog','Display','off','ConstraintTolerance',constr_tol);
N = size(l,3);
Ain_coll = []; 
bin_coll = []; 
success = 0;
outbound = 0;
coll = 0;

for k = 1: k_hor
    [violation,min_dist,viol_constr] = CheckCollSoftDMPC(prev_p(:,k),l,n,k,E1,rmin,order);
    if (any(violation))
        N_violation = sum(viol_constr);
        if (k==1 && min_dist < rmin - 0.05) %already violated constraint
            p = [];
            v = [];
            a = [];
            coll = 1;
            return;
            
        elseif (k==1)
            continue;
        end     
        [Ainr,binr,prev_dist] = CollConstrSoftDMPC2(prev_p(:,k),po,vo,n,k,l,rmin,A,A_initp,E1,E2,order,viol_constr);
        Ain_coll = [Ainr diag(prev_dist)];
%         Ain_coll = [Ainr diag(prev_dist)];
        bin_coll = binr;
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

if (any(violation)) % In case of collisions, we relax the constraint with slack variable
    % Add dimensions for slack variable
    Q = [Q zeros(3*K,N_violation);
         zeros(N_violation,3*K) zeros(N_violation,N_violation)];
    R = [R zeros(3*K,N_violation);
         zeros(N_violation,3*K) zeros(N_violation,N_violation)];
    S = [S zeros(3*K,N_violation);
         zeros(N_violation,3*K) zeros(N_violation,N_violation)];
    A = [A zeros(3*K,N_violation);
         zeros(N_violation,3*K) zeros(N_violation,N_violation)];
    Delta = [Delta zeros(3*K,N_violation);
         zeros(N_violation,3*K) zeros(N_violation,N_violation)];
    bin_total = [bin_coll; repmat((pmax)',K,1) - A_initp*([po';vo']); zeros(N_violation,1); repmat(-(pmin)',K,1) + A_initp*([po';vo']); zeros(N_violation,1)];
    ao_1 = [ao zeros(1,3*(K-1)+N_violation)];
    A_initp_aug = [A_initp; zeros(N_violation,6)];
    
    % add bound on the relaxation variable
    ub = [ub; zeros(N_violation,1)];
    lb = [lb; -0.05*ones(N_violation,1)];

    % Linear penalty on collision constraint relaxation
    f_eps = term*[zeros(3*K,1); ones(N_violation,1)]';
    
    % Quadratic penalty on collision constraint relaxation
    EPS = 1*10^0*[zeros(3*K,3*K) zeros(3*K,N_violation);
           zeros(N_violation,3*K) eye(N_violation,N_violation)];
       
    f = -2*([repmat((pf)',K,1); zeros(N_violation,1)]'*Q*A - (A_initp_aug*([po';vo']))'*Q*A + ao_1*S*Delta) + f_eps ;
    
else % case of no collisions, we don't even add collision constraints
    bin_total = [bin_coll; repmat((pmax)',K,1) - A_initp*([po';vo']); repmat(-(pmin)',K,1) + A_initp*([po';vo'])];
    ao_1 = [ao zeros(1,3*(K-1))];
    f = -2*(repmat((pf)',K,1)'*Q*A - (A_initp*([po';vo']))'*Q*A + ao_1*S*Delta) ;
    EPS = zeros(3*K,3*K); 
end

Ain_total = [Ain_coll; A; -A];
H = 2*(A'*Q*A+ Delta'*S*Delta + R + EPS);
tries = 0;
x = [];
%Solve and propagate states
while(~success && tries < 30)
    [x,fval,exitflag] = quadprog(H,f',Ain_total,bin_total,Aeq,beq,lb,ub,[],options);
    if (exitflag == -6)
        % Weird non-convex flag may appear, even though the problem is
        % very well defined as a convex problem
        % fix: increase constraint tolerance and retry solving
        p = [];
        v = [];
        a = [];
        fprintf("Exitflag was -6 in Bound \n")
        constr_tol = 2*constr_tol;
        options.ConstraintTolerance = constr_tol;
        success = 0;
        tries = tries + 1;
        continue  
    elseif (~isempty(x))
        % everything was good, return the solution
        a = x(1:3*K);
        [p,v] = propStatedmpc(po,vo,a,A_initp, A_p, A_v);
        p = vec2mat(p,3)';
        v = vec2mat(v,3)';
        a = vec2mat(a,3)';
        success = 1;
        if (~is_inbounds(p(:,1),pmin,pmax))
            success = 0;
            outbound = 1;
        end 
%         if violation % extract the value of the slack variable (not used atm)
%          epsilon = x(3*K+1:end);
% %     fprintf("min epsilon = %.4f e-3 \n",1000*min(epsilon))
%         end
        return
        
    elseif isempty(x)
        p = [];
        v = [];
        a = [];
        success = 0;
        if (~any(violation))
            fprintf("Couldn't solve in no violation case, retrying with higher tolerance \n");
            constr_tol = 2*constr_tol;
            options.ConstraintTolerance = constr_tol;
            tries = tries + 1;
            continue
        end
        fprintf("Retrying with more relaxed bound \n");
        lb(3*K+1:end) = 2*lb(3*K+1:end);
        term = term*2;
        f_eps = term*[zeros(3*K,1); ones(N_violation,1)]';
        f = -2*([repmat((pf)',K,1); zeros(N_violation,1)]'*Q*A - (A_initp_aug*([po';vo']))'*Q*A + ao_1*S*Delta) + f_eps ;
        tries = tries + 1;
        continue
    end
end
if ~success
    hola = 1;
end

end