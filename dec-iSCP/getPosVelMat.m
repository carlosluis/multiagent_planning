function Aaug = getPosVelMat(h,K)

% Kinematic model A,b matrices
A = [1 0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

b = [h^2/2*eye(3);
     h*eye(3)];
 
prev_row = zeros(6,3*K); % For the first iteration of constructing matrix Ain

% Build matrix to convert acceleration to position
for k = 1:K
    add_b = [zeros(size(b,1),size(b,2)*(k-1)) b zeros(size(b,1),size(b,2)*(K-k))];
    new_row = A*prev_row + add_b;   
    prev_row = new_row; 
end

Aaug = [new_row; [zeros(3,3*(K-1)) eye(3)]; [eye(3) zeros(3,3*(K-1))]];

end