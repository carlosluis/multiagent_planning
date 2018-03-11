function [p,v,a] = singleiSCP(po, pf, h, K, pmin,pmax, l)

prev_p = initSolution(po,pf,h,K);
epsilon = 2; %to be tuned
tol = 0;
H = eye(3*K);
a_lim = 1; %Maximum acc of 1 m/s^2
ub = a_lim*ones(3*K,1);
lb = -ub; 
i = 1;

r_min = 0.5; %minimum radius in meters

% Kinematic model A,b matrices
A = [1 0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

b = [h^2/2*eye(3);
     h*eye(3)];
 
bin = []; % local inequality constraint for any given drone
Ain = []; % local inequality constraint variable
diff_mat = []; % 
prev_row = zeros(6,3*K); % For the first iteration of constructing matrix Ain

 
% Build matrix to convert acceleration to position
for k = 1:K
    add_b = [zeros(size(b,1),size(b,2)*(k-1)) b zeros(size(b,1),size(b,2)*(K-k))];
    new_row = A*prev_row + add_b;   
    Ain = [Ain; new_row(1:3,:)];
    prev_row = new_row; 
end



addConstr = [];

while i <= 11
    newConstrCount = 0; 
    Ain_total = [];
    bin_total = [];
    for k = 1:K
        violation = CheckforColl (prev_p(:,k),l,k,r_min);
        if (ismember(k,addConstr))
            [Ainr, binr] = CollConstr(prev_p(:,k),po,k, l, Ain,r_min);
            Ain_total = [Ain_total; Ainr];
            bin_total = [bin_total; binr];
        
        elseif (newConstrCount==0 && violation)
            [Ainr, binr] = CollConstr(prev_p(:,k-1),po,k-1, l, Ain,r_min);
            Ain_total = [Ain_total; Ainr];
            bin_total = [bin_total; binr];  
            addConstr = [addConstr k];
            newConstrCount = newConstrCount + 1;
        end           
    end
    
    % Setup the QP
    
    Ain_total = [Ain_total; Ain; -Ain];
    bin_total = [bin_total; repmat(pmax',K,1); repmat(-pmin',K,1)];
    
    [Aeq, beq] = EqConstr(pf-po,K,h);
    
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