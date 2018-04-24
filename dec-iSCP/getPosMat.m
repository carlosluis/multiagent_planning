function Apos = getPosMat(h,K)

% Kinematic model A,b matrices
A = [1 0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

b = [h^2/2*eye(3);
     h*eye(3)];
 
Apos = zeros(3*K,3*K); % local inequality constraint variable
prev_row = zeros(6,3*K); % For the first iteration of constructing matrix Ain
idx=1;
% Build matrix to convert acceleration to position
for k = 1:K
    add_b = [zeros(size(b,1),size(b,2)*(k-1)) b zeros(size(b,1),size(b,2)*(K-k))];
    new_row = A*prev_row + add_b;   
    Apos(idx:idx+2,:) = new_row(1:3,:);
    prev_row = new_row; 
    idx = idx+3;
end