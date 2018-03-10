function [Ain_total, bin_total] = IneqConstr(p, l, h, pmin, pmax)

r_min = 0.9; %minimum radius in meters

% Kinematic model A,b matrices
A = [1 0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

b = [h^2/2*eye(3);
     h*eye(3)];
 
K = size(l,2); % number of time steps of the trajectory K = T/h +1
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

% Add the boundary constraints pmin and pmax
Ain_total = [Ain; -Ain];
bin_total = [repmat(pmax',K,1); repmat(-pmin',K,1)];

if (~isempty(l))
    for i = 1:size(l,3) %Iterate through the number of obstacles (other agents)
        pj = l(:,:,i); %position vector of the i-th neighbour over k iterations

        for k = 1:K % Iterate through all time steps of the trajectory
            dist = norm(p(:,k)-pj(:,k)); %distance at time step k
            diff = (p-pj)'; % Transpose of the difference

            % Right side of inequality constraint (bin)
            r = dist*(r_min - dist + (p(:,k) - pj(:,k))'*p(:,k));
            bin = [bin; r];

            % Construct diagonal matrix with vector difference
            diff_mat = [diff_mat; zeros(1,3*(k-1)) diff(1,:) zeros(1,3*(K-k))];
        end

        % Update the ineq. constraints matrix and vector
        Ain_total =  [Ain_total; -diff_mat*Ain];
        bin_total = [bin_total; -bin];

        % Reset vectors to calculate next obstacle
        bin = [];
        diff_mat = [];  
    end
end
end