function tol = maxDeviation(p, prev_p)

K = length(p)/3;

for k = 1:K
    dist(k) = norm(p(:,k) - prev_p(:,k)); %distance between iters at time step k 
end

tol = max(dist);

end