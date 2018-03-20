function [Ain_total, bin_total] = AllCollConstr(p,K,rmin,l,A)
Ain_total = [];
bin_total = [];
bin = [];
diff_mat = [];

if (~isempty(l))
    for i = 1:size(l,3) %Iterate through the number of obstacles (other agents)
        pj = l(:,:,i); %position vector of the i-th neighbour over k iterations

        for k = 1:K % Iterate through all time steps of the trajectory
            dist = norm(p(:,k)-pj(:,k)); %distance at time step k
            diff = (p(:,k)-pj(:,k))'; % Transpose of the difference

            % Right side of inequality constraint (bin)
            r = dist*(rmin - dist + (p(:,k) - pj(:,k))'*p(:,k)/dist) - (p(:,k) - pj(:,k))'*p(:,1);
            bin = [bin; r];

            % Construct diagonal matrix with vector difference
            diff_mat = [diff_mat; zeros(1,3*(k-1)) diff zeros(1,3*(K-k))];
        end

        % Update the ineq. constraints matrix and vector
        Ain_total =  [Ain_total; -diff_mat*A];
        bin_total = [bin_total; -bin];

        % Reset vectors to calculate next obstacle
        bin = [];
        diff_mat = [];  
    end
end
end