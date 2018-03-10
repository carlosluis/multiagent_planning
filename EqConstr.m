function [Aeq, beq] = EqConstr(pf,K,h)

% Kinematic model A,b matrices
A = [1 0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

b = [h^2/2*eye(3);
     h*eye(3)];
 
beq = [pf' ; zeros(3,1); zeros(3,1)];
prev_row = zeros(6,3*K);
Aeq = prev_row;

for k = 1:K % Iterate through all time steps of the trajectory
    % Left side of the equality constraint (Aeq)
    add_b = [zeros(size(b,1),size(b,2)*(k-1)) b zeros(size(b,1),size(b,2)*(K-k))];
    new_row = A*prev_row + add_b;   
    prev_row = new_row;
end

Aeq = [new_row; [zeros(3,3*(K-1)) eye(3)]];

end
