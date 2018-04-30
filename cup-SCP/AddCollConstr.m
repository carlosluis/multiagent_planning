function [Ain_total, bin_total] = AddCollConstr(p,po,K,rmin,A,E1,E2,order)

N = size(p,3);
Ain_total = zeros(K*N*(N-1)/2,3*K*N);
bin_total = zeros(K*N*(N-1)/2,1);
l = 1;

for i = 1:N-1 % i-th vehicle
    pi = p(:,:,i);
    for j = i+1:N % j-th neighbour
        pj = p(:,:,j);
        for k = 1:K % k-th time step of trajectory
            dist = norm(E1*(pi(:,k)-pj(:,k)),order);
            diff = (E2*(pi(:,k)-pj(:,k)).^(order-1));
            
            % Right side of inequality constraint (bin)
            r = dist^(order-1)*(rmin - dist) + diff'*diff - diff'*(po(:,:,i)'-po(:,:,j)');
            bin = r;
            
            % Construct diagonal matrix with vector difference
            diff_mat = [zeros(1,3*K*(i-1)) zeros(1,3*(k-1))  diff' zeros(1,3*(K-k)) zeros(1,3*K*(j-i-1)) zeros(1,3*(k-1)) -diff' zeros(1,3*(K-k)) zeros(1,3*K*(N-j))];
            
            % Update the ineq. constraints matrix and vector
            Ain_total(l,:) =  -diff_mat*A;
            bin_total(l) = -bin;  
            l = l+1;    
        end
    end   
end
        
    