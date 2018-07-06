function [violation,min_dist,viol_constr] = CheckCollSoftDMPC (p,l,n,k,E1,rmin,order)
idx = 1;
N = size(l,3);
violation = zeros(N,1);
viol_constr = zeros(N,1);
for i = 1:N %Iterate through the number of obstacles (other agents)
    if(i~=n)
        pj = l(:,:,i); %position vector of the i-th neighbour over k iterations
        dist(idx) = norm(E1*(p-pj(:,k)),order);
        violation(i) = (dist(idx) < rmin);
        viol_constr(i) = dist(idx) < rmin*(1+k/15);
        idx = idx + 1;
    end  
end
min_dist = min(dist);
end
