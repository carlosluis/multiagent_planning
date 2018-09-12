function Delta = getDeltaMat(k_hor)
K = k_hor - 1;
Delta = [eye(3) zeros(3,3*(K))];
b = [-eye(3) eye(3)];

for k = 1:K
    new_row = [zeros(3,3*(k-1)) b zeros(3,size(b,2)*(K-k)/2)];
    Delta = [Delta; new_row];
end
    