function [Ain_total, bin_total] = AddCollConstr(p,po,K,rmin,A)

Ain_total = [];
bin_total = [];
bin = [];
diff_mat = [];
N = size(p,3);

for i = 1:N-1 % i-th vehicle
    pi = p(:,:,i);
    for j = i+1:N % j-th neighbour
        pj = p(:,:,j);
        
        for k = 1:K % k-th time step of trajectory
            dist = norm(pi(:,k)-pj(:,k));
            diff = pi(:,k)-pj(:,k);
            eta = diff./dist;
            
            % Right side of inequality constraint (bin)
            r = rmin - dist + eta'*diff - eta'*(po(:,:,i)'-po(:,:,j)');
            bin = [bin;r];
            
            % Construct diagonal matrix with vector difference
            diff_mat = [diff_mat; zeros(1,3*K*(i-1)) zeros(1,3*(k-1))  eta' zeros(1,3*(K-k)) zeros(1,3*K*(j-i-1)) zeros(1,3*(k-1)) -eta' zeros(1,3*(K-k)) zeros(1,3*K*(N-j))];
            
            % Update the ineq. constraints matrix and vector
            Ain_total =  [Ain_total; -diff_mat*A];
            bin_total = [bin_total; -bin];  
            
            % Reset vectors to calculate next obstacle
            bin = [];
            diff_mat = [];        
        end 
    end   
end
        
    