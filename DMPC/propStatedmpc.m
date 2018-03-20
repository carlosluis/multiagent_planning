function [p,v] = propStatedmpc(po, vo,  a, h)

% Kinematic model A,b matrices
A = [1 0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

b = [h^2/2*eye(3);
     h*eye(3)];

K = length(a)/3;

prev_row = zeros(6,3*K); % For the first iteration of constructing matrix Ain
A_p = [];
A_v = [];
A_initp = [];
A_init = eye(6);

% Build matrix to convert acceleration to position
for k = 1:(K)
    add_b = [zeros(size(b,1),size(b,2)*(k-1)) b zeros(size(b,1),size(b,2)*(K-k))];
    new_row = A*prev_row + add_b;   
    A_p = [A_p; new_row(1:3,:)];
    A_v = [A_v; new_row(4:6,:)];
    prev_row = new_row;
    A_init = A*A_init;
    A_initp = [A_initp; A_init(1:3,:)];
end

new_p = A_p*a + repmat(po',K,1);
new_v = A_v*a;

p = new_p;
v = new_v;
end