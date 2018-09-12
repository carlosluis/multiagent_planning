function [p,v,a,success,outbound,coll] = solveSoftDMPCbound_debug(po,pf,vo,ao,n,h,l,K,rmin,pmin,pmax,alim,A,A_initp,A_p,A_v,Delta,Q1,S1,E1,E2,order,term)

k_hor = size(l,2);
ub = alim*ones(3*K,1);
lb = -ub; 
prev_p = l(:,:,n);
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
        end   
        [Ainr,binr,prev_dist] = CollConstrSoftDMPC(prev_p(:,k),po,vo,n,k,l,rmin,A,A_initp,E1,E2,order,viol_constr);
        Ain_coll = [Ainr diag(prev_dist)];
        bin_coll = binr;
        break;
    end       
end

% Plot the separating planes resulting from the distance constraint
% linearization

colors = distinguishable_colors(N);

if (any(violation))
   viol_idx = find(violation==1);
   num_coll = length(viol_idx);
   figure(1)
%    set(gcf, 'Position', get(0, 'Screensize'));
   [az,el] = view;
   clf
   for i=1:num_coll
       p1 = prev_p(:,k);
       p2 = l(:,k,viol_idx(i));
       diff = (p2 - p1);
       w = null(diff');
       new_p = l(:,k,viol_idx(i)) - E1^(-1)*diff/norm(diff)*rmin;
       [P,Q] = meshgrid(-1:1); % Provide a gridwork (you choose the size)
       X = new_p(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
       Y = new_p(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
       Z = new_p(3)+w(3,1)*P+w(3,2)*Q;
       surf(X,Y,Z)
       hold on;
       xlim([-3,3])
       ylim([-3,3])
       zlim([-3,3])
%        plot3(l(1,:,viol_idx(i)),l(2,:,viol_idx(i)),l(3,:,viol_idx(i)),...
%            'o','Color',colors(viol_idx(i),:),'Linewidth',2)
       plot3(prev_p(1,:),prev_p(2,:),prev_p(3,:),...
           'o','Color',colors(n,:),'Linewidth',1)
       plot3(p1(1),p1(2),p1(3),'o','Color',colors(n,:),'Linewidth',4);
       plot3(p2(1),p2(2),p2(3),'o','Color',colors(viol_idx(i),:),'Linewidth',4);
       plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)])
       hold off
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
        if (any(violation))
            figure(1)
            hold on;
            plot3(p(1,:),p(2,:),p(3,:),...
            '*','Color',colors(n,:),'Linewidth',2)
            plot3(p(1,k),p(2,k),p(3,k),...
               '*','Color',colors(n,:),'Linewidth',8)
            plot3(pf(1),pf(2),pf(3),'s','Color',colors(n,:),'Linewidth',4,'markers',10)
            view([az,el]);
            hold off
            pause(.01)
        end
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