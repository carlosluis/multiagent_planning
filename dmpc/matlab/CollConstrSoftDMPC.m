function [Ain_total, bin_total,prev_dist] = CollConstrSoftDMPC(p,po,vo,n, k, l,rmin, Ain,A_initp,E1,E2,order,violation)
N_obs = size(l,3);
N_violation = sum(violation);
Ain_total = zeros(N_violation,3*size(l,2));
bin_total = zeros(N_violation,1);
prev_dist = zeros(N_violation,1);
idx = 1;
k_ctr = k;

if (~isempty(l))
    for i = 1:N_obs %Iterate through the number of obstacles (other agents)
        if(i~=n && violation(i))
            pj = l(:,:,i); %poisition vector of the i-th neighbour over k iterations
            K = size(pj,2);

            dist = norm(E1*(p-pj(:,k)),order); %distance at time step k
            diff = (E2*(p-pj(:,k)).^(order-1))'; % Transpose of the difference
            
            prev_dist(idx) = dist^(order-1);
            % Right side of inequality constraint (bin)
            r = dist^(order-1)*((rmin) - dist + diff*p/(dist^(order-1))) - diff*A_initp(3*(k_ctr-1)+1:3*k_ctr,:)*[po';vo'];

            % Construct diagonal matrix with vector difference
            diff_mat = [zeros(1,3*(k_ctr-1)) diff zeros(1,3*(K-k_ctr))];

            % Update the ineq. constraints matrix and vector
            Ain_total(idx,:) =  -diff_mat*Ain;
            bin_total(idx,:) = -r; 
            idx = idx + 1;
        end
    end
end