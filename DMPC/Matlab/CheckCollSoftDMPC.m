function [violation,min_dist] = CheckCollSoftDMPC (p,l,n,k,E1,rmin,order)
violation = false;
idx = 1;
for i = 1:size(l,3) %Iterate through the number of obstacles (other agents)
    if(i~=n)
        pj = l(:,:,i); %position vector of the i-th neighbour over k iterations
        dist(idx) = norm(E1*(p-pj(:,k)),order);
        violation = (violation) || (dist(idx) < rmin);
        idx = idx + 1;
    end  
end
min_dist = min(dist);
end
