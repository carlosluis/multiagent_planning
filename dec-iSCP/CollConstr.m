function [Ain_total, bin_total] = CollConstr(p,po, k, l,Ain,rmin,E1,E2,order)
N_obs = size(l,3);
Ain_total = zeros(N_obs,size(Ain,1));
bin_total = zeros(N_obs,1);
if (~isempty(l))
    for i = 1:size(l,3) %Iterate through the number of obstacles (other agents)
        pj = l(:,:,i); %position vector of the i-th neighbour over k iterations
        K = size(pj,2);
        
        dist = norm(E1*(p-pj(:,k)),order); %distance at time step k
        diff = (E2*(p-pj(:,k)).^(order-1))'; % Transpose of the difference

        % Right side of inequality constraint (bin)
        r = dist^(order-1)*(rmin - dist + diff*p/(dist^(order-1))) - diff*po';

        % Construct diagonal matrix with vector difference
        diff_mat = [zeros(1,3*(k-2)) diff zeros(1,3*(K-(k-1)))];

        % Update the ineq. constraints matrix and vector
        Ain_total(i,:) =  -diff_mat*Ain;
        bin_total(i) = -r;
    end
end